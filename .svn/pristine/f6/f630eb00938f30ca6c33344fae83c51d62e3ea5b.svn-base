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


#include "Hamiltonian/ParticleOnCylinderInterlayerPfaffian.h"
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


#define M1_12 0.08333333333333333


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderInterlayerPfaffian::ParticleOnCylinderInterlayerPfaffian(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum,
										   double ratio,AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();

  //S2 = 1/2(Nu + Nd) + sum (m1,m2) c_m1u^+ c_m2d^+ c_m2u c_m1d
  this->OneBodyUpUp = new Complex[this->NbrLzValue];
  this->OneBodyDownDown = new Complex[this->NbrLzValue];
  Complex OneHalf(0.5, 0.0);
  for (int i = 0; i <= this->MaxMomentum; i++)
	this->OneBodyUpUp[i] = OneHalf;
  for (int i = 0; i <= this->MaxMomentum; i++)
	this->OneBodyDownDown[i] = OneHalf;

  this->OneBodyUpDown = new Complex[this->NbrLzValue];

  //set Infinity to large value to get the psi sector
  //by not allowing more than 1 particles in orbitals 0 and 2S
  Complex Infinity(0.0, 0.0);
  Complex Zero(0.0, 0.0);
  this->OneBodyUpDown[0] = Infinity;
  this->OneBodyUpDown[this->MaxMomentum] = Infinity;

  for (int i = 1; i < this->MaxMomentum; i++)
	this->OneBodyUpDown[i] = Zero;
    
  this->EnergyShift = 0.0;

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
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnCylinderInterlayerPfaffian::~ParticleOnCylinderInterlayerPfaffian() 
{
  delete[] this->InteractionFactors12;
  delete[] this->M1Value12;
  delete[] this->M2Value12;
  delete[] this->M3Value12;
  delete[] this->M4Value12;
  delete[] this->M5Value12;
  delete[] this->M6Value12;
  delete[] this->InteractionFactors32;
  delete[] this->M1Value32;
  delete[] this->M2Value32;
  delete[] this->M3Value32;
  delete[] this->M4Value32;
  delete[] this->M5Value32;
  delete[] this->M6Value32;
  delete[] this->InteractionFactorsIntra;
  delete[] this->M1ValueIntra;
  delete[] this->M2ValueIntra;
  delete[] this->M3ValueIntra;
  delete[] this->M4ValueIntra;
  delete[] this->InteractionFactorsInter;
  delete[] this->M1ValueInter;
  delete[] this->M2ValueInter;
  delete[] this->M3ValueInter;
  delete[] this->M4ValueInter;

  delete[] this->OneBodyUpUp;
  delete[] this->OneBodyUpDown;
  delete[] this->OneBodyDownDown;


  if (this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      this->FastMultiplicationFlag = false;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnCylinderInterlayerPfaffian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnCylinderInterlayerPfaffian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionFactors()
{
  int Pos;
  int m4;
  double MaxCoefficient;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {

      //2-body intra
     
      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 < m3))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientIntra(m1, m2, m3, m4);
		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsIntra = 0;
      this->M1ValueIntra = new int [Pos];
      this->M2ValueIntra = new int [Pos];
      this->M3ValueIntra = new int [Pos];
      this->M4ValueIntra = new int [Pos];

      this->InteractionFactorsIntra = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 < m3))
		     {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsIntra[this->NbrInteractionFactorsIntra] = TmpCoefficient[Pos];
		           this->M1ValueIntra[this->NbrInteractionFactorsIntra] = m1;
		           this->M2ValueIntra[this->NbrInteractionFactorsIntra] = m2;
		           this->M3ValueIntra[this->NbrInteractionFactorsIntra] = m3;
		           this->M4ValueIntra[this->NbrInteractionFactorsIntra] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsIntra;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction 2-body, S=1 = " << this->NbrInteractionFactorsIntra << endl;
     cout << "====================================" << endl;


      //2-body inter
     
      Pos = 0;    
      MaxCoefficient = 0.0;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientInter(m1, m2, m3, m4);
		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsInter = 0;
      this->M1ValueInter = new int [Pos];
      this->M2ValueInter = new int [Pos];
      this->M3ValueInter = new int [Pos];
      this->M4ValueInter = new int [Pos];

      this->InteractionFactorsInter = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		     {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsInter[this->NbrInteractionFactorsInter] = TmpCoefficient[Pos];
		           this->M1ValueInter[this->NbrInteractionFactorsInter] = m1;
		           this->M2ValueInter[this->NbrInteractionFactorsInter] = m2;
		           this->M3ValueInter[this->NbrInteractionFactorsInter] = m3;
		           this->M4ValueInter[this->NbrInteractionFactorsInter] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsInter;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction 2-body, S=0 = " << this->NbrInteractionFactorsInter << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;




//----------------------------------- 3BODY ----------------------------------------------------
      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient12 = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient12[Pos] = this->EvaluateInteractionCoefficient12(m1, m2, m3, m4, m5, m6);
		       if (MaxCoefficient < Norm(TmpCoefficient12[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient12[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactors12 = 0;
      this->M1Value12 = new int [Pos];
      this->M2Value12 = new int [Pos];
      this->M3Value12 = new int [Pos];
      this->M4Value12 = new int [Pos];
      this->M5Value12 = new int [Pos];
      this->M6Value12 = new int [Pos];

      this->InteractionFactors12 = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		     {
		       if (Norm(TmpCoefficient12[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactors12[this->NbrInteractionFactors12] = TmpCoefficient12[Pos];
		           this->M1Value12[this->NbrInteractionFactors12] = m1;
		           this->M2Value12[this->NbrInteractionFactors12] = m2;
		           this->M3Value12[this->NbrInteractionFactors12] = m3;
		           this->M4Value12[this->NbrInteractionFactors12] = m4;
		           this->M5Value12[this->NbrInteractionFactors12] = m5;
		           this->M6Value12[this->NbrInteractionFactors12] = m6;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors12;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction S=1/2 = " << this->NbrInteractionFactors12 << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient12;

//----------------------------------------------------------------------------------------------------------

      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient32 = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < m2; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 < m5))
  	             {
 		       TmpCoefficient32[Pos] = this->EvaluateInteractionCoefficient32(m1, m2, m3, m4, m5, m6);
		       if (MaxCoefficient < Norm(TmpCoefficient32[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient32[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactors32 = 0;
      this->M1Value32 = new int [Pos];
      this->M2Value32 = new int [Pos];
      this->M3Value32 = new int [Pos];
      this->M4Value32 = new int [Pos];
      this->M5Value32 = new int [Pos];
      this->M6Value32 = new int [Pos];

      this->InteractionFactors32 = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < m2; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 < m5))
		     {
		       if (Norm(TmpCoefficient32[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactors32[this->NbrInteractionFactors32] = TmpCoefficient32[Pos];
		           this->M1Value32[this->NbrInteractionFactors32] = m1;
		           this->M2Value32[this->NbrInteractionFactors32] = m2;
		           this->M3Value32[this->NbrInteractionFactors32] = m3;
		           this->M4Value32[this->NbrInteractionFactors32] = m4;
		           this->M5Value32[this->NbrInteractionFactors32] = m5;
		           this->M6Value32[this->NbrInteractionFactors32] = m6;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors32;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction S=3/2 = " << this->NbrInteractionFactors32 << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient32;



    }
  else //bosons
    {
      cout << "Bosons require more testing. " << endl;
      exit(2);

      //2-body intra
     
      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= m3))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientIntra(m1, m2, m3, m4) * this->NumberOfPermutations12(m1, m2) * this->NumberOfPermutations12(m3, m4);
		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsIntra = 0;
      this->M1ValueIntra = new int [Pos];
      this->M2ValueIntra = new int [Pos];
      this->M3ValueIntra = new int [Pos];
      this->M4ValueIntra = new int [Pos];

      this->InteractionFactorsIntra = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= m3))
		     {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsIntra[this->NbrInteractionFactorsIntra] = TmpCoefficient[Pos];
		           this->M1ValueIntra[this->NbrInteractionFactorsIntra] = m1;
		           this->M2ValueIntra[this->NbrInteractionFactorsIntra] = m2;
		           this->M3ValueIntra[this->NbrInteractionFactorsIntra] = m3;
		           this->M4ValueIntra[this->NbrInteractionFactorsIntra] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsIntra;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction 2-body, S=1 = " << this->NbrInteractionFactorsIntra << endl;
     cout << "====================================" << endl;


      //2-body inter
     
      Pos = 0;    
      MaxCoefficient = 0.0;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientInter(m1, m2, m3, m4);
		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsInter = 0;
      this->M1ValueInter = new int [Pos];
      this->M2ValueInter = new int [Pos];
      this->M3ValueInter = new int [Pos];
      this->M4ValueInter = new int [Pos];

      this->InteractionFactorsInter = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
 	       {
	         m4 = m1 + m2 - m3;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		     {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsInter[this->NbrInteractionFactorsInter] = TmpCoefficient[Pos];
		           this->M1ValueInter[this->NbrInteractionFactorsInter] = m1;
		           this->M2ValueInter[this->NbrInteractionFactorsInter] = m2;
		           this->M3ValueInter[this->NbrInteractionFactorsInter] = m3;
		           this->M4ValueInter[this->NbrInteractionFactorsInter] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsInter;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction 2-body, S=0 = " << this->NbrInteractionFactorsInter << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;


//--------------------------------- 3BODY ----------------------------------------------

      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient12 = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 <= m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient12[Pos] = this->EvaluateInteractionCoefficientBosons12(m1, m2, m3, m4, m5, m6) * this->NumberOfPermutations12(m1, m2) * this->NumberOfPermutations12(m5, m6);
		       if (MaxCoefficient < Norm(TmpCoefficient12[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient12[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactors12 = 0;
      this->M1Value12 = new int [Pos];
      this->M2Value12 = new int [Pos];
      this->M3Value12 = new int [Pos];
      this->M4Value12 = new int [Pos];
      this->M5Value12 = new int [Pos];
      this->M6Value12 = new int [Pos];

      this->InteractionFactors12 = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 <= m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		     {
		       if (Norm(TmpCoefficient12[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactors12[this->NbrInteractionFactors12] = TmpCoefficient12[Pos];
		           this->M1Value12[this->NbrInteractionFactors12] = m1;
		           this->M2Value12[this->NbrInteractionFactors12] = m2;
		           this->M3Value12[this->NbrInteractionFactors12] = m3;
		           this->M4Value12[this->NbrInteractionFactors12] = m4;
		           this->M5Value12[this->NbrInteractionFactors12] = m5;
		           this->M6Value12[this->NbrInteractionFactors12] = m6;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors12;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction S=1/2 = " << this->NbrInteractionFactors12 << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient12;

//----------------------------------------------------------------------------------------------------------

      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient32 = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= m2; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 <= m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= m5))
  	             {
 		       TmpCoefficient32[Pos] = this->EvaluateInteractionCoefficientBosons32(m1, m2, m3, m4, m5, m6) * this->NumberOfPermutations32(m1, m2, m3) * this->NumberOfPermutations32(m4, m5, m6);
		       if (MaxCoefficient < Norm(TmpCoefficient32[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient32[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactors32 = 0;
      this->M1Value32 = new int [Pos];
      this->M2Value32 = new int [Pos];
      this->M3Value32 = new int [Pos];
      this->M4Value32 = new int [Pos];
      this->M5Value32 = new int [Pos];
      this->M6Value32 = new int [Pos];

      this->InteractionFactors32 = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= m2; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 <= m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= m5))
		     {
		       if (Norm(TmpCoefficient32[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactors32[this->NbrInteractionFactors32] = TmpCoefficient32[Pos];
		           this->M1Value32[this->NbrInteractionFactors32] = m1;
		           this->M2Value32[this->NbrInteractionFactors32] = m2;
		           this->M3Value32[this->NbrInteractionFactors32] = m3;
		           this->M4Value32[this->NbrInteractionFactors32] = m4;
		           this->M5Value32[this->NbrInteractionFactors32] = m5;
		           this->M6Value32[this->NbrInteractionFactors32] = m6;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors32;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction S=3/2 = " << this->NbrInteractionFactors32 << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient32;

    }
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficientIntra(int m1, int m2, int m3, int m4)
{
  Complex Zero(0.0, 0.0);
  return Zero; 
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficientInter(int m1, int m2, int m3, int m4)
{

  Complex Coefficient(0.0, 0.0);
  
  if (m1==m4)
    Coefficient.Re += 1.0;

  return (Coefficient);
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// m5 = fifth index
// m6 = sixth index
// return value = numerical coefficient

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficient12(int m1, int m2, int m3, int m4, int m5, int m6)
{

  //Complex Zero(0,0);
  //return Zero;

  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m6 - m4 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

     GaussianExp = exp(-Xr * Xr - Xs * Xs - Xr * Xs) * exp(-Xrp * Xrp - Xsp * Xsp - Xrp * Xsp);


     //Coefficient.Re = GaussianExp * (Xr - Xs) * (Xrp - Xsp) * (1.0 +  (2.0 * Xr + Xs) * (2.0 * Xs + Xr) * (2.0 * Xrp + Xsp) * (2.0 * Xsp + Xrp) + (Xr + Xs) * (Xrp + Xsp) );

     Coefficient.Re = GaussianExp * ( 0.5 * (Xr - Xs) * (Xrp - Xsp) + 1.5 * (Xr * Xr - Xs * Xs) * (Xrp * Xrp - Xsp * Xsp) ); 

     return (-Coefficient * 16.0 * sqrt(M_PI) * sqrt(3.0 * M_PI) /(2.0 * M_PI * this->Ratio * this->NbrLzValue));
 
}

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficient32(int m1, int m2, int m3, int m4, int m5, int m6)
{
  //Complex Zero(0,0);
  //return Zero; 

  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m6 - m4 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

     GaussianExp = Xr * Xr + Xs * Xs + Xr * Xs;
     Coefficient.Re = exp(-GaussianExp) * (Xr - Xs) * (2.0 * Xr + Xs) * (2.0 * Xs + Xr);
     Coefficient.Im = 0.0;

     GaussianExp = Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;
     Coefficient.Re *= (exp(-GaussianExp) * (Xrp - Xsp) * (2.0 * Xrp + Xsp) * (2.0 * Xsp + Xrp));
     Coefficient.Im = 0.0;
     return (-Coefficient * 16.0 * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
 
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// m5 = fifth index
// m6 = sixth index
// return value = numerical coefficient

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficientBosons12(int m1, int m2, int m3, int m4, int m5, int m6)
{
  Complex Coefficient(0,0);

  return Coefficient; 
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// m5 = fifth index
// m6 = sixth index
// return value = numerical coefficient

Complex ParticleOnCylinderInterlayerPfaffian::EvaluateInteractionCoefficientBosons32(int m1, int m2, int m3, int m4, int m5, int m6)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m6 - m4 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

  GaussianExp = Xr * Xr + Xs * Xs + Xr * Xs;
  Coefficient.Re = exp(-GaussianExp);
  Coefficient.Im = 0.0;

  GaussianExp = Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;
  Coefficient.Re *= exp(-GaussianExp);
  Coefficient.Im = 0.0;
  
  return (Coefficient * (2.0/3.0) * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
 
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnCylinderInterlayerPfaffian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = this->EnergyShift;
  double Coefficient, Coefficient2;

  if (this->FastMultiplicationFlag == false)
    {
      //cout<<"Flagfalse"<<endl;

      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      int m5;
      int m6;
      int TmpA[3];
      int TmpC[3];
      int TmpSpinA[3];
      int TmpSpinC[3];
      int TmpA2b[2];
      int TmpC2b[2];
      int TmpSpinA2b[2];
      int TmpSpinC2b[2];


      Complex TmpInteraction;
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();


      for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	{

          //UPUP

	  TmpC2b[0] = this->M1ValueIntra[j];
	  TmpC2b[1] = this->M2ValueIntra[j];
	  TmpA2b[0] = this->M3ValueIntra[j];
	  TmpA2b[1] = this->M4ValueIntra[j];

	  TmpSpinC2b[0] = 0;
	  TmpSpinC2b[1] = 0;
	  TmpSpinA2b[0] = 0;
	  TmpSpinA2b[1] = 0;

	  TmpInteraction = this->InteractionFactorsIntra[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }


          //DOWNDOWN

	  TmpC2b[0] = this->M1ValueIntra[j];
	  TmpC2b[1] = this->M2ValueIntra[j];
	  TmpA2b[0] = this->M3ValueIntra[j];
	  TmpA2b[1] = this->M4ValueIntra[j];

	  TmpSpinC2b[0] = 1;
	  TmpSpinC2b[1] = 1;
	  TmpSpinA2b[0] = 1;
	  TmpSpinA2b[1] = 1;

	  TmpInteraction = this->InteractionFactorsIntra[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }


	} //j

      for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	{

          //UPDOWN

	  TmpC2b[0] = this->M1ValueInter[j];
	  TmpC2b[1] = this->M2ValueInter[j];
	  TmpA2b[0] = this->M3ValueInter[j];
	  TmpA2b[1] = this->M4ValueInter[j];

	  TmpSpinC2b[0] = 0;
	  TmpSpinC2b[1] = 1;
	  TmpSpinA2b[0] = 0;
	  TmpSpinA2b[1] = 1;

	  TmpInteraction = this->InteractionFactorsInter[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }

	} //j


      for (int j = 0; j < this->NbrInteractionFactors12; ++j) 
	{

          //UPUPDOWN 

	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
	  TmpSpinA[2] = 0;
          
	  TmpInteraction = this->InteractionFactors12[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }


          //DOWNDOWNUP 

	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
	  TmpSpinA[2] = 1;
          
	  TmpInteraction = this->InteractionFactors12[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<"DDU "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }

	} //j

      for (int j = 0; j < this->NbrInteractionFactors32; ++j) 
	{

          //UPUPUP

	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;
	  TmpSpinA[2] = 0;
          
	  TmpInteraction = this->InteractionFactors32[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }


          //DOWNDOWNDOWN

	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;
	  TmpSpinA[2] = 1;
          
	  TmpInteraction = this->InteractionFactors32[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<"DDU "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }

	} //j

      for (int i = firstComponent; i < LastComponent; ++i)
	{
	 vDestination[i] += this->EnergyShift * vSource[i];
	}

      //One-body terms
      double TmpUp, TmpDown; 
      for (int i = firstComponent; i < LastComponent; ++i)
	  { 
            Complex TmpDiagonal(0, 0);
 	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      {
                TmpUp = TmpParticles->AduAu(i, j); 
                TmpDown = TmpParticles->AddAd(i, j); 
                if (this->OneBodyUpUp != 0)
			TmpDiagonal += this->OneBodyUpUp[j] * TmpUp;
                if (this->OneBodyDownDown != 0)
			TmpDiagonal += this->OneBodyDownDown[j] * TmpDown;
                if (this->OneBodyUpDown != 0)
			TmpDiagonal += this->OneBodyUpDown[j] * TmpUp * TmpDown;
	      }

           //TmpUp = TmpParticles->AduAu(i, this->MaxMomentum); 
           //TmpDown = TmpParticles->AddAd(i, this->MaxMomentum); 
           //if ((TmpUp==0) && (TmpDown==0))
           //    TmpDiagonal += 1000.0;
           //if ((TmpUp==1) && (TmpDown==1))
           //   TmpDiagonal += 1000.0;
  
           //TmpUp = TmpParticles->AduAu(i, 0); 
           //TmpDown = TmpParticles->AddAd(i, 0); 
           //if ((TmpUp==0) && (TmpDown==0))
           //    TmpDiagonal += 1000.0;
           //if ((TmpUp==1) && (TmpDown==1))
           //    TmpDiagonal += 1000.0;

 
	    vDestination[i] += (TmpDiagonal * vSource[i]);
	  }


/*
         for (int i = firstComponent; i < LastComponent; ++i)
	  { 
 	    for (int j = 0; j < this->MaxMomentum; ++j) 
	      {
                Index = TmpParticles->AduAd(i, j, j+1, Coefficient); 
	        if (Index < Dim)
 		  vDestination[Index] += Coefficient * 100.0 * vSource[i];

                Index = TmpParticles->AddAu(i, j+1, j, Coefficient); 
	        if (Index < Dim)
 		  vDestination[Index] += Coefficient * 100.0 * vSource[i];
	      }
           }

*/


      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
          ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (j = 0; j < TmpNbrInteraction; ++j)
               { 
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
                  if (TmpIndexArray[j] < i)
                     vDestination[i] += vSource[TmpIndexArray[j]] * Conj(TmpCoefficientArray[j]);
               }
	      vDestination[i] += Shift * vSource[i];
	    }


        //One-body terms
        double TmpUp, TmpDown; 
        for (int i = firstComponent; i < LastComponent; ++i)
	  { 
            Complex TmpDiagonal(0, 0);;
 	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      {
                TmpUp = TmpParticles->AduAu(i, j); 
                TmpDown = TmpParticles->AddAd(i, j); 
                if (this->OneBodyUpUp != 0)
			TmpDiagonal += this->OneBodyUpUp[j] * TmpUp;
                if (this->OneBodyDownDown != 0)
			TmpDiagonal += this->OneBodyDownDown[j] * TmpDown;
                if (this->OneBodyUpDown != 0)
			TmpDiagonal += this->OneBodyUpDown[j] * TmpUp * TmpDown;
	      }
	    vDestination[i] += (TmpDiagonal * vSource[i]);
	  }

        
          delete TmpParticles;
	}
      else
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      for (j = 0; j < TmpNbrInteraction; ++j)
               {
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
                  if (TmpIndexArray[j] < i)
                     vDestination[i] += vSource[TmpIndexArray[j]] * Conj(TmpCoefficientArray[j]);
               }
	      vDestination[i] += Shift * vSource[i];
	      ++Pos;
	    }


        //One-body terms
        double TmpUp, TmpDown; 
        for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	  { 
            Complex TmpDiagonal(0, 0);;
 	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      {
                TmpUp = TmpParticles->AduAu(i, j); 
                TmpDown = TmpParticles->AddAd(i, j); 
                if (this->OneBodyUpUp != 0)
			TmpDiagonal += this->OneBodyUpUp[j] * TmpUp;
                if (this->OneBodyDownDown != 0)
			TmpDiagonal += this->OneBodyDownDown[j] * TmpDown;
                if (this->OneBodyUpDown != 0)
			TmpDiagonal += this->OneBodyUpDown[j] * TmpUp * TmpDown;
	      }
	    vDestination[i] += (TmpDiagonal * vSource[i]);
	  }


	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
          int m5;
          int m6;
          int TmpA[3];
          int TmpC[3];
          int TmpSpinA[3];
          int TmpSpinC[3];
          int TmpA2b[2];
          int TmpC2b[2];
          int TmpSpinA2b[2];
          int TmpSpinC2b[2];


	  Complex TmpInteraction;

	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		


    		  for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
			{

		        //UPUP

			TmpC2b[0] = this->M1ValueIntra[j];
	 		TmpC2b[1] = this->M2ValueIntra[j];
	  		TmpA2b[0] = this->M3ValueIntra[j];
	  		TmpA2b[1] = this->M4ValueIntra[j];

			TmpSpinC2b[0] = 0;
	  		TmpSpinC2b[1] = 0;
	  		TmpSpinA2b[0] = 0;
	  		TmpSpinA2b[1] = 0;

	  		TmpInteraction = this->InteractionFactorsIntra[j];
	  		for (int i = firstComponent; i < LastComponent; ++i)
	   		 {
	      		Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
              		if (Coefficient2 != 0.0)
                		{
                   		Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
  	    		       if (Index <= i)
            		        {
				      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
            		          //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
            		          //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
            		          if (Index < i)
            		            vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
            		        }
            		    } 
	    		}


         		 //DOWNDOWN

			  TmpC2b[0] = this->M1ValueIntra[j];
	 		 TmpC2b[1] = this->M2ValueIntra[j];
	  		TmpA2b[0] = this->M3ValueIntra[j];
	  		TmpA2b[1] = this->M4ValueIntra[j];

			  TmpSpinC2b[0] = 1;
	  		TmpSpinC2b[1] = 1;
	  		TmpSpinA2b[0] = 1;
	  		TmpSpinA2b[1] = 1;

			  TmpInteraction = this->InteractionFactorsIntra[j];
	 		 for (int i = firstComponent; i < LastComponent; ++i)
	    		{
	     		 Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
              		if (Coefficient2 != 0.0)
                		{
                   		Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
		  	           if (Index <= i)
                		    {
				      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                		      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                		      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                		      if (Index < i)
                		        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                		    }
                		} 
	    		}


			} //j

		      for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
			{

		          //UPDOWN

			  TmpC2b[0] = this->M1ValueInter[j];
	 		 TmpC2b[1] = this->M2ValueInter[j];
	  		TmpA2b[0] = this->M3ValueInter[j];
	  		TmpA2b[1] = this->M4ValueInter[j];

			  TmpSpinC2b[0] = 0;
	  		TmpSpinC2b[1] = 1;
	  		TmpSpinA2b[0] = 0;
	  		TmpSpinA2b[1] = 1;

			  TmpInteraction = this->InteractionFactorsInter[j];
	 		 for (int i = firstComponent; i < LastComponent; ++i)
	    		{
	      		Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
		              if (Coefficient2 != 0.0)
                		{
                		   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
  	        		   if (Index <= i)
                		    {
				      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                		      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                		      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                		      if (Index < i)
                		        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                		    }
                		} 
	    		}

			} //j




	      for (int j = 0; j < this->NbrInteractionFactors12; ++j) 
		{

	          //UPUPDOWN 

		  TmpC[0] = this->M1Value12[j];
		  TmpC[1] = this->M2Value12[j];
		  TmpC[2] = this->M3Value12[j];
		  TmpA[0] = this->M4Value12[j];
		  TmpA[1] = this->M5Value12[j];
		  TmpA[2] = this->M6Value12[j];

		  TmpSpinC[0] = 0;
		  TmpSpinC[1] = 0;
		  TmpSpinC[2] = 1;
		  TmpSpinA[0] = 1;
		  TmpSpinA[1] = 0;
		  TmpSpinA[2] = 0;
          
		  TmpInteraction = this->InteractionFactors12[j];
		  for (int i = firstComponent; i < LastComponent; ++i)
		    {
		      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
	              if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
	  	           if (Index <= i)
	                    {
			      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
	                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
	                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
	                      if (Index < i)
	                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
	                    }
	                } 
		    }


	          //DOWNDOWNUP 

		  TmpC[0] = this->M1Value12[j];
		  TmpC[1] = this->M2Value12[j];
		  TmpC[2] = this->M3Value12[j];
		  TmpA[0] = this->M4Value12[j];
		  TmpA[1] = this->M5Value12[j];
		  TmpA[2] = this->M6Value12[j];

		  TmpSpinC[0] = 1;
		  TmpSpinC[1] = 1;
		  TmpSpinC[2] = 0;
		  TmpSpinA[0] = 0;
		  TmpSpinA[1] = 1;
		  TmpSpinA[2] = 1;
	          
		  TmpInteraction = this->InteractionFactors12[j];
		  for (int i = firstComponent; i < LastComponent; ++i)
		    {
		      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
	              if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
	  	           if (Index <= i)
	                    {
			      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
	                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
	                      //cout<<"DDU "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
	                      if (Index < i)
	                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
	                    }
	                } 
		    }
	
		} //j

	      for (int j = 0; j < this->NbrInteractionFactors32; ++j) 
		{

	          //UPUPUP

		  TmpC[0] = this->M1Value32[j];
		  TmpC[1] = this->M2Value32[j];
		  TmpC[2] = this->M3Value32[j];
		  TmpA[0] = this->M4Value32[j];
		  TmpA[1] = this->M5Value32[j];
		  TmpA[2] = this->M6Value32[j];

		  TmpSpinC[0] = 0;
		  TmpSpinC[1] = 0;
		  TmpSpinC[2] = 0;
		  TmpSpinA[0] = 0;
		  TmpSpinA[1] = 0;
		  TmpSpinA[2] = 0;
          
		  TmpInteraction = this->InteractionFactors32[j];
		  for (int i = firstComponent; i < LastComponent; ++i)
		    {
		      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
	              if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
	  	           if (Index <= i)
	                    {
			      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
	                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
	                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
	                      if (Index < i)
	                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
	                    }
	                } 
		    }


	          //DOWNDOWNDOWN

		  TmpC[0] = this->M1Value32[j];
		  TmpC[1] = this->M2Value32[j];
		  TmpC[2] = this->M3Value32[j];
		  TmpA[0] = this->M4Value32[j];
		  TmpA[1] = this->M5Value32[j];
		  TmpA[2] = this->M6Value32[j];

		  TmpSpinC[0] = 1;
		  TmpSpinC[1] = 1;
		  TmpSpinC[2] = 1;
		  TmpSpinA[0] = 1;
		  TmpSpinA[1] = 1;
		  TmpSpinA[2] = 1;
          
		  TmpInteraction = this->InteractionFactors32[j];
		  for (int i = firstComponent; i < LastComponent; ++i)
		    {
		      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
	              if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
	  	           if (Index <= i)
	                    {
			      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
	                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
	                      //cout<<"DDU "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
	                      if (Index < i)
	                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
	                    }
	                } 
		    }
	
		} //j

		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    vDestination[i] += this->EnergyShift * vSource[i];
		  }

              //One-body terms
              double TmpUp, TmpDown; 
              for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	        { 
                  Complex TmpDiagonal(0, 0);;
 	          for (int j = 0; j <= this->MaxMomentum; ++j) 
	            {
                       TmpUp = TmpParticles->AduAu(i, j); 
                       TmpDown = TmpParticles->AddAd(i, j); 
                       if (this->OneBodyUpUp != 0)
	  		   TmpDiagonal += this->OneBodyUpUp[j] * TmpUp;
                       if (this->OneBodyDownDown != 0)
			   TmpDiagonal += this->OneBodyDownDown[j] * TmpDown;
                       if (this->OneBodyUpDown != 0)
			   TmpDiagonal += this->OneBodyUpDown[j] * TmpUp * TmpDown;
	            }
	          vDestination[i] += (TmpDiagonal * vSource[i]);
	        }

                
	      }
	  delete TmpParticles;
	}
   }


  return vDestination;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnCylinderInterlayerPfaffian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  int m5;
  int m6;
  int TmpA[3];
  int TmpC[3];
  int TmpSpinA[3];
  int TmpSpinC[3];
  int TmpA2b[2];
  int TmpC2b[2];
  int TmpSpinA2b[2];
  int TmpSpinC2b[2];

  int LastComponent = lastComponent + firstComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  cout<<"Partial fast memory permanent "<<endl;

  for (int i = firstComponent; i < LastComponent; ++i)
    {

	for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	  {

	          //UPUP

		  TmpC2b[0] = this->M1ValueIntra[j];
		  TmpC2b[1] = this->M2ValueIntra[j];
		  TmpA2b[0] = this->M3ValueIntra[j];
		  TmpA2b[1] = this->M4ValueIntra[j];

		  TmpSpinC2b[0] = 0;
		  TmpSpinC2b[1] = 0;
		  TmpSpinA2b[0] = 0;
		  TmpSpinA2b[1] = 0;

                  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
                              ++Memory;
	                      ++this->NbrInteractionPerComponent[i];                    
	                    }
	                } 
		    

	          //DOWNDOWN

		  TmpC2b[0] = this->M1ValueIntra[j];
		  TmpC2b[1] = this->M2ValueIntra[j];
		  TmpA2b[0] = this->M3ValueIntra[j];
		  TmpA2b[1] = this->M4ValueIntra[j];

		  TmpSpinC2b[0] = 1;
		  TmpSpinC2b[1] = 1;
		  TmpSpinA2b[0] = 1;
		  TmpSpinA2b[1] = 1;

                  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
                              ++Memory;
	                      ++this->NbrInteractionPerComponent[i];                    
	                    }
	                } 
		
	
	} //j


	for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	  {

	          //UPUP

		  TmpC2b[0] = this->M1ValueInter[j];
		  TmpC2b[1] = this->M2ValueInter[j];
		  TmpA2b[0] = this->M3ValueInter[j];
		  TmpA2b[1] = this->M4ValueInter[j];

		  TmpSpinC2b[0] = 0;
		  TmpSpinC2b[1] = 1;
		  TmpSpinA2b[0] = 0;
		  TmpSpinA2b[1] = 1;

                  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
                              ++Memory;
	                      ++this->NbrInteractionPerComponent[i];                    
	                    }
	                } 
		    	
	} //j


      for (int j = 0; j < this->NbrInteractionFactors12; ++j) 
	{

	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	}  //j

      for (int j = 0; j < this->NbrInteractionFactors32; ++j) 
	{

	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	}  //j  
    }
  Memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * TmpParticles->GetHilbertSpaceDimension() + 
	    Memory *  (sizeof (int) + sizeof(Complex)));
  delete TmpParticles;
  return Memory;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void ParticleOnCylinderInterlayerPfaffian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  int NbrTranslation;
  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();  
  int LastComponent = firstComponent + lastComponent;
  int m1, m2, m3, m4, m5, m6;
  int TmpC[3];
  int TmpA[3]; 
  int TmpSpinC[3];
  int TmpSpinA[3]; 
  int TmpC2b[2];
  int TmpA2b[2]; 
  int TmpSpinC2b[2];
  int TmpSpinA2b[2]; 

  int SumIndices;
  int Pos;
  int ReducedNbrInteractionFactors;
  int PosIndex = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++PosIndex;
      PosMod = this->FastMultiplicationStep - PosMod;
    }

  cout<<"Partial enable fast permanent "<<endl;

  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)  
    {

      this->InteractionPerComponentIndex[PosIndex] = new int [this->NbrInteractionPerComponent[PosIndex]];
      this->InteractionPerComponentCoefficient[PosIndex] = new Complex [this->NbrInteractionPerComponent[PosIndex]];      
      TmpIndexArray = this->InteractionPerComponentIndex[PosIndex];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[PosIndex];
      Pos = 0;

	for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	  {

	          //UPUP

		  TmpC2b[0] = this->M1ValueIntra[j];
		  TmpC2b[1] = this->M2ValueIntra[j];
		  TmpA2b[0] = this->M3ValueIntra[j];
		  TmpA2b[1] = this->M4ValueIntra[j];

		  TmpSpinC2b[0] = 0;
		  TmpSpinC2b[1] = 0;
		  TmpSpinA2b[0] = 0;
		  TmpSpinA2b[1] = 0;

		  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
			       TmpIndexArray[Pos] = Index;
			       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsIntra[j];
			       ++Pos;
	                    }
	                } 
		 

	          //DOWNDOWN

		  TmpC2b[0] = this->M1ValueIntra[j];
		  TmpC2b[1] = this->M2ValueIntra[j];
		  TmpA2b[0] = this->M3ValueIntra[j];
		  TmpA2b[1] = this->M4ValueIntra[j];

		  TmpSpinC2b[0] = 1;
		  TmpSpinC2b[1] = 1;
		  TmpSpinA2b[0] = 1;
		  TmpSpinA2b[1] = 1;

		  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
	       			TmpIndexArray[Pos] = Index;
	       			TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsIntra[j];
	       			++Pos;
	                    }
	                } 
		 
	} //j


	for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	  {

	          //UPUP

		  TmpC2b[0] = this->M1ValueInter[j];
		  TmpC2b[1] = this->M2ValueInter[j];
		  TmpA2b[0] = this->M3ValueInter[j];
		  TmpA2b[1] = this->M4ValueInter[j];

		  TmpSpinC2b[0] = 0;
		  TmpSpinC2b[1] = 1;
		  TmpSpinA2b[0] = 0;
		  TmpSpinA2b[1] = 1;

		  Coefficient2 = TmpParticles->ProdA(i, TmpA2b, TmpSpinA2b, 2);
	          if (Coefficient2 != 0.0)
	                {
	                   Index = TmpParticles->ProdAd(TmpC2b, TmpSpinC2b, 2, Coefficient);
	  	           if (Index <= i)
	                    {
			       TmpIndexArray[Pos] = Index;
			       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsInter[j];
			       ++Pos;
	                    }
	                } 
		 
	} //j

      for (int j = 0; j < this->NbrInteractionFactors12; ++j) 
	{
	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors12[j];
	       ++Pos;
              }
           }


	  TmpC[0] = this->M1Value12[j];
	  TmpC[1] = this->M2Value12[j];
	  TmpC[2] = this->M3Value12[j];
	  TmpA[0] = this->M4Value12[j];
	  TmpA[1] = this->M5Value12[j];
	  TmpA[2] = this->M6Value12[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors12[j];
	       ++Pos;
              }
           }

 
	} //j

      for (int j = 0; j < this->NbrInteractionFactors32; ++j) 
	{
	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors32[j];
	       ++Pos;
              }
           }


	  TmpC[0] = this->M1Value32[j];
	  TmpC[1] = this->M2Value32[j];
	  TmpC[2] = this->M3Value32[j];
	  TmpA[0] = this->M4Value32[j];
	  TmpA[1] = this->M5Value32[j];
	  TmpA[2] = this->M6Value32[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors32[j];
	       ++Pos;
              }
           }

 
	} //j

     ++PosIndex;
   }

  delete TmpParticles;
}

// Get the number of permutations of annihilation/creation indices c_n1 c_n2 c_n3 for bosons

int ParticleOnCylinderInterlayerPfaffian:: NumberOfPermutations32(int n1, int n2, int n3)
{
  if ((n1 != n2) && (n2 != n3) && (n3 != n1) ) 
    return 6;
  else if ( ((n1==n2) && (n1 != n3)) || ((n1==n3) && (n1 != n2)) || ((n2==n3) && (n2 != n1)) ) 
    return 3;
  else
    return 1;
}

// Get the number of permutations of annihilation/creation indices c_n1 c_n2 c_n3 for bosons

int ParticleOnCylinderInterlayerPfaffian:: NumberOfPermutations12(int n1, int n2)
{
  if (n1 != n2)
    return 2;
  else
    return 1;
}
