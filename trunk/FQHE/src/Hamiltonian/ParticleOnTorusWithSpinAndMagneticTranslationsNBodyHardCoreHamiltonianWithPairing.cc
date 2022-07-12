////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//         n-body hardcore interaction and magnetic translations              //
//                                                                            //
//                        last modification : 23/10/2015                      //
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


#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include "GeneralTools/ArrayTools.h"

#include "Architecture/AbstractArchitecture.h"

#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333

// default constructor
//

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing::ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing()
{
  this->PairingAmplitude = 0.0;
}


// constructor from pseudopotentials
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrNBody = value of the n (i.e. the n-body interaction)
// intraUpSpinNBodyInteraction = strength of the n-body interaction for particles with spin up
// intraDownSpinNBodyInteraction = strength of the n-body interaction for particles with spin down
// interSpinNBodyInteraction = strength of the n-body interaction for particles with different spincomponent
// nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
// pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
// nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
// pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
// nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
// pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
// pairing =  amplitude of the pariring term
// spinFluxUp = additional inserted flux for spin up
// spinFluxDown = additional inserted flux for spin down
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing::ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody, 
																	       double intraUpSpinNBodyInteraction, 
																	       double intraDownSpinNBodyInteraction, 
																	       double interSpinNBodyInteraction, 
																	       int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
																	       int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
																	       int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
																	       double pairing, double spinFluxUp, double spinFluxDown, 
																	       AbstractArchitecture* architecture, long memory, char* precalculationFileName, double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;//true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  


  this->NBodyValue = nbrNBody;
  this->IntraUpSpinNBodyInteractionStrength = intraUpSpinNBodyInteraction;
  this->IntraDownSpinNBodyInteractionStrength = intraDownSpinNBodyInteraction;
  this->InterSpinNBodyInteractionStrength = interSpinNBodyInteraction;

  this->PairingAmplitude = pairing;
  this->MaxNBody = this->NBodyValue;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NbrSpinSectors = new int [this->MaxNBody + 1];
  for (int i = 0; i < this->MaxNBody; ++i)
    this->NBodyFlags[i] = false;
  this->NBodyFlags[this->MaxNBody] = true;

  this->FullTwoBodyFlag = false;
  this->NbrPseudopotentialsUpUp = nbrPseudopotentialsUpUp;
  if (pseudopotentialsUpUp != 0)
    {
      this->PseudopotentialsUpUp = new double[this->NbrPseudopotentialsUpUp];
      for (int i = 0; i < this->NbrPseudopotentialsUpUp; ++i)
	this->PseudopotentialsUpUp[i] = pseudopotentialsUpUp[i];
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->PseudopotentialsUpUp = 0;
    }
  this->NbrPseudopotentialsDownDown = nbrPseudopotentialsDownDown;
  if (pseudopotentialsDownDown != 0)
    {
      this->PseudopotentialsDownDown = new double[this->NbrPseudopotentialsDownDown];
      for (int i = 0; i < this->NbrPseudopotentialsDownDown; ++i)
	this->PseudopotentialsDownDown[i] = pseudopotentialsDownDown[i];
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->PseudopotentialsDownDown = 0;
    }
  this->NbrPseudopotentialsUpDown = nbrPseudopotentialsUpDown;
  if (pseudopotentialsUpDown != 0)
    {
      this->PseudopotentialsUpDown = new double[this->NbrPseudopotentialsUpDown];
      for (int i = 0; i < this->NbrPseudopotentialsUpDown; ++i)
	this->PseudopotentialsUpDown[i] = pseudopotentialsUpDown[i];
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->PseudopotentialsUpDown = 0;
    }

  this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpUp;
  if (this->NbrPseudopotentialsDownDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsDownDown;
  if (this->NbrPseudopotentialsUpDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpDown;
  if (this->MaxNbrPseudopotentials > 0)
    {
      this->LaguerrePolynomials = new Polynomial[this->MaxNbrPseudopotentials];
      for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
	this->LaguerrePolynomials[i] = LaguerrePolynomial(i);
    }
  else
    {
      this->LaguerrePolynomials = 0;
    }

  this->OneBodyInteractionFactorsupup = 0;
  if(oneBodyPotentielUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsupup[i] = oneBodyPotentielUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if(oneBodyPotentielDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsdowndown[i] = oneBodyPotentielDownDown[i];
    } 
  this->OneBodyInteractionFactorsupdown = 0;
  if(oneBodyPotentielUpDown != 0)
    {
      this->OneBodyInteractionFactorsupdown = new Complex[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsupdown[i] = oneBodyPotentielUpDown[i];
	}
    } 

  this->EvaluateExponentialFactors();

  this->NbrEntryPrecalculatedInteractionCoefficients1 = this->MaxMomentum;
  for (int i = 1; i < this->NBodyValue; ++i)
    this->NbrEntryPrecalculatedInteractionCoefficients1 *= this->MaxMomentum;
  this->PrecalculatedInteractionCoefficients = new double* [this->NbrEntryPrecalculatedInteractionCoefficients1];
  this->NbrEntryPrecalculatedInteractionCoefficients2 = this->NBodyValue*(2*this->NBodyValue - 1);
  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
    {      
      this->PrecalculatedInteractionCoefficients[m1] = new double [this->NbrEntryPrecalculatedInteractionCoefficients2];
    }
  char* InteractionCoefficientFileName = new char [512];
  sprintf (InteractionCoefficientFileName, "%dbodydelta_interactioncoefficient_2s_%d_ratio_%.10f.dat", this->NBodyValue, this->NbrLzValue, this->Ratio);
  if (IsFile(InteractionCoefficientFileName))
    {
      ifstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "cannot open " << InteractionCoefficientFileName << endl;
	}
      else
	{
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    ReadBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	  File.close();
	}
    }
  else
    {
      ofstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::out);
      if (!File.is_open())
	{
	  cout << "cannot create " << InteractionCoefficientFileName << endl;
	}
      else
	{  
	    
	  int** TmpIndices = new int* [this->NbrEntryPrecalculatedInteractionCoefficients1];
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    {
	      TmpIndices[m1] = new int [this->NBodyValue];
	      int Tmp = m1;
	      for (int i = 0; i < this->NBodyValue; ++i)
		{
		  TmpIndices[m1][i] = Tmp % this->MaxMomentum;
		  Tmp /= this->MaxMomentum;
		}
	    }
	  	  
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    {
	      int TmpIndex = 0;
	      for (int j = 0; j < 2*this->NBodyValue - 1; ++j)
		{
		  int momentumTransfer = j - this->NBodyValue + 1;
		  double* Coefficient = this->EvaluateInteractionCoefficientCreation(TmpIndices[m1], momentumTransfer);
		  for (int g = 0; g < this->NBodyValue; ++g)
		    {
		      this->PrecalculatedInteractionCoefficients[m1][TmpIndex] = Coefficient[g];
		      ++TmpIndex;
		    }
		  delete[] Coefficient;
		}
	      WriteBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	    }

	  
	  delete[] InteractionCoefficientFileName;
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    delete[] TmpIndices[m1];
	  delete[] TmpIndices;
	  File.close();
	}
    }
  

  this->EvaluateInteractionFactors();

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024l)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1l << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1l << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  cout << endl;
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

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing::~ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing() 
{
}

// evaluate all interaction factors
//   

void ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;

  if (this->InterSpinNBodyInteractionStrength == 0.0)
    {
      this->GetIndices(true);
    }
  else
    {
      this->GetIndices(false);
    }

  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout << "ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonianWithPairing::EvaluateInteractionFactors not implemented for  fermions" << endl;
    }
  else
    {      
      this->NBodyInteractionFactors = new Complex***[this->MaxNBody + 1];
      for (int TmpNBody = 0; TmpNBody <= this->MaxNBody; ++TmpNBody)
	{
	  if (this->NBodyFlags[TmpNBody] == true)
	    {
// 	      for (int j = 0; j <  this->NbrSpinSectors[TmpNBody]; ++j)
// 		{
// 		  int TmpNbrSumSector = this->NbrNBodySpinMomentumSectorSum[TmpNBody][j];
// 		  int** LinearizedNBodySectorIndicesPerSum = new int* [TmpNbrSumSector];
// 		  for (int i = 0; i < TmpNbrSumSector; ++i)
// 		    {
// 		      LinearizedNBodySectorIndicesPerSum[i] = new int [this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i]];
// 		      for (int j1 = 0; j1 < this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i]; ++j1)
// 			{
// 			  int Tmp = 0;
// 			  for (int k = TmpNBody - 1; k >= 0; --k)
// 			    {
// 			      Tmp *= this->MaxMomentum;
// 			      Tmp += this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i][(j1 * this->NBodyValue) + k];
// 			    }
// 			  LinearizedNBodySectorIndicesPerSum[i][j1] = Tmp;
// 			}
// 		      SortArrayUpOrdering<int>(LinearizedNBodySectorIndicesPerSum[i], this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i]);
// 		      for (int j1 = 0; j1 < this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i]; ++j1)
// 			{
// 			  int Tmp = LinearizedNBodySectorIndicesPerSum[i][j1];
// 			  for (int k = 0; k < TmpNBody; ++k)
// 			    { 
// 			      this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][j][i][(j1 * this->NBodyValue) + k] = Tmp % this->MaxMomentum;
// 			      Tmp /= this->MaxMomentum;
// 			    }	  
// 			}
// 		    }
// 		}

	      int TmpMaxSpinSector = this->NbrSpinSectors[TmpNBody] - 1;
	      this->NBodyInteractionFactors[TmpNBody] = new Complex**[this->NbrSpinSectors[TmpNBody]];
	      int TmpNbrSumSector = this->NbrNBodySpinMomentumSectorSum[TmpNBody][0];
	      this->NBodyInteractionFactors[TmpNBody][0] = new Complex*[TmpNbrSumSector];
	      this->NBodyInteractionFactors[TmpNBody][TmpMaxSpinSector] = new Complex*[TmpNbrSumSector];
	      for (int j = 0; j < TmpNbrSumSector; ++j)
		{
		  this->NBodyInteractionFactors[TmpNBody][0][j] = new Complex[this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][0][j] 
									      * this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][0][j]];
		  this->NBodyInteractionFactors[TmpNBody][TmpMaxSpinSector][j] = new Complex[this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][TmpMaxSpinSector][j]
											 * this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][TmpMaxSpinSector][j]];
		  int Index = 0;  	  
		  int Lim = this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][0][j];
		  for (int j1 = 0; j1 < Lim; ++j1)
		    {
		      int* TmpNIndices = &(this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][0][j][j1 * TmpNBody]);
		      for (int j2 = 0; j2 < Lim; ++j2)
			{
			  int* TmpMIndices = &(this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][0][j][j2 * TmpNBody]);
			  double Tmp = this->EvaluateInteractionNIndexSymmetrizedCoefficient(TmpNIndices, TmpMIndices);
			  this->NBodyInteractionFactors[TmpNBody][0][j][Index] = Tmp * this->IntraUpSpinNBodyInteractionStrength;
			  this->NBodyInteractionFactors[TmpNBody][TmpMaxSpinSector][j][Index] = Tmp * this->IntraDownSpinNBodyInteractionStrength;
			  TotalNbrInteractionFactors += 2;
			  ++Index;		  		  
			}
		    }
 		}
	      for (int i = 1; i < TmpMaxSpinSector; ++i)
		{
		  int TmpNbrSumSector = this->NbrNBodySpinMomentumSectorSum[TmpNBody][i];
		  this->NBodyInteractionFactors[TmpNBody][i] = new Complex*[TmpNbrSumSector];
		  for (int j = 0; j < TmpNbrSumSector; ++j)
		    {
		      this->NBodyInteractionFactors[TmpNBody][i][j] = new Complex[this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][i][j] 
										  * this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][i][j]];
		      int Index = 0;  	  
		      int Lim = this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][i][j];
		      for (int j1 = 0; j1 < Lim; ++j1)
			{
			  int* TmpNIndices = &(this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][i][j][j1 * TmpNBody]);
			  for (int j2 = 0; j2 < Lim; ++j2)
			    {
			      int* TmpMIndices = &(this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][i][j][j2 * TmpNBody]);
			      double Tmp = this->EvaluateInteractionNIndexTwoSetSymmetrizedCoefficient(TmpNIndices, TmpMIndices, TmpNBody - i, i);
			      this->NBodyInteractionFactors[TmpNBody][i][j][Index] = ((double) TmpNBody) * Tmp * this->InterSpinNBodyInteractionStrength;
			      TotalNbrInteractionFactors++;
			      ++Index;		  		  
			    }
			}
		    }
		}
	    }
	}
    }
  if (this->FullTwoBodyFlag == true)
    {
      double MaxCoefficient = 0.0;
      
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  // upup-upup 
	  this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsupup[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  // downdown-downdown
	  this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
	  MaxCoefficient = 0.0;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
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
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										      this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										      this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsdowndown[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  // updown-updown
	  this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
	  MaxCoefficient = 0.0;
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
						 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsupdown[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  cout << "warning, pairing is not implemented for fermions" << endl;
	}
      else
	{
	  // upup-upup 
	  this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsupup[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  // downdown-downdown
	  this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
	  MaxCoefficient = 0.0;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
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
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										      this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										      this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
											this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsdowndown[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  // upup-downdown
	  this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
	  MaxCoefficient = 0.0;
	  double* FakePseudopotentials = new double [1];
	  FakePseudopotentials[0] = this->PairingAmplitude;
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
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 1, FakePseudopotentials, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 1, FakePseudopotentials, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 1, FakePseudopotentials,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 1, FakePseudopotentials,
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 1, FakePseudopotentials, 
										      this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 1, FakePseudopotentials, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 1, FakePseudopotentials, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 1, FakePseudopotentials, 
											this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown));
		      if (m1 == m2)		    
			TmpCoefficient *= 0.5;
		      if (m3 == m4)		    
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsupupdowndown[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsupupdowndown[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	  delete[] FakePseudopotentials;
	  // updown-updown
	  this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
	  MaxCoefficient = 0.0;
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										      this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										      this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
						 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
											this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
			  ++TotalNbrNonZeroInteractionFactors;
			}
		      else
			{
			  this->InteractionFactorsupdown[i][Index] = 0.0;
			}
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	}
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "nbr non-zero interaction = " << TotalNbrNonZeroInteractionFactors << endl;
  cout << "====================================" << endl;
}

