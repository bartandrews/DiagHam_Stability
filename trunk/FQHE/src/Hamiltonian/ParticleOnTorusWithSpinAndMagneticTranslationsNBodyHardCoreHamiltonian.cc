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


#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian.h"
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

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian()
{
  this->Particles = 0;
  this->LzMax = 0;
  this->NbrLzValue = 0;
  this->NbrParticles = 0;
  this->MomentumModulo = 0;
  this->Ratio = 0.0;
  this->InvRatio = 0.0;
  this->MaxMomentum = 0;
  this->XMomentum = 0;
  
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;//true;
  this->Ratio = 0;  
  this->InvRatio = 0;
  this->SpinFluxUp = 0;
  this->SpinFluxDown = 0;
  this->HamiltonianShift = 0.0;
  this->Architecture = 0;
  this->PrecalculationShift = 0;  

  this->NbrPseudopotentialsUpUp = 0;
  this->PseudopotentialsUpUp = 0;
  this->NbrPseudopotentialsDownDown = 0;
  this->NbrPseudopotentialsUpDown = 0;
  this->PseudopotentialsUpDown = 0;
  
  this->MaxNbrPseudopotentials = 0;
  this->LaguerrePolynomials = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;

  this->NBodyValue = 0;

  this->IntraUpSpinNBodyInteractionStrength = 0.0;
  this->IntraDownSpinNBodyInteractionStrength = 0.0;
  this->InterSpinNBodyInteractionStrength = 0.0;

  this->FullTwoBodyFlag = false;
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
// spinFluxUp = additional inserted flux for spin up
// spinFluxDown = additional inserted flux for spin down
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody, 
																	       double intraUpSpinNBodyInteraction, 
																	       double intraDownSpinNBodyInteraction, 
																	       double interSpinNBodyInteraction, 
																	       int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
																	       int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
																	       int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
																	       double spinFluxUp, double spinFluxDown, 
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

ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::~ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian() 
{
  if (this->PseudopotentialsUpUp != 0)
    delete[] this->PseudopotentialsUpUp;
  if (this->PseudopotentialsDownDown != 0)
    delete[] this->PseudopotentialsDownDown;
  if (this->PseudopotentialsUpDown != 0)
    delete[] this->PseudopotentialsUpDown;
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
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
      cout << "ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateInteractionFactors not implemented for  fermions" << endl;
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

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// spinFluxM1 = additional inserted flux for m1
// spinFluxM2 = additional inserted flux for m2
// spinFluxM3 = additional inserted flux for m3
// spinFluxM4 = additional inserted flux for m4
// return value = numerical coefficient

double ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													      double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - (((double) (m1-m3)) + spinFluxM1 - spinFluxM3) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = ((double) (m1 - m4)) + spinFluxM1 - spinFluxM4;
  double N1;
  double Q2;
  double Precision;
  double TmpInteraction;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0* PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
 	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue) + spinFluxM1 - spinFluxM4;
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i=0; i< nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  //Normalize per flux (gives correct energy scale for 2-particle problem)
  return (Sum / ((double) this->NbrLzValue));
}

// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation or the annihilation operators only) for each integer modulo the NBodyValue
//
// mIndices = array that contains the creation indices
// tmpSum = sum of all creation indices
// momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
// return value = array of numerical coefficients  

double* ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateInteractionCoefficientCreation(int* mIndices, int momentumTransfer)
{
  double* Coefficient = new double [this->NBodyValue];
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  int tmpSum = 0;
  for (int i = 0; i < this->NBodyValue; ++i)
    tmpSum += mIndices[i];
  
  int* momFactor = new int [this->NBodyValue - 1];
  int* TmpIndices = new int [(this->NBodyValue - 1)];
  int* countIter = new int[(this->NBodyValue - 1)];
  for (int i = 0; i < this->NBodyValue - 1; ++i)
    momFactor[i] = this->NBodyValue * mIndices[i] -  tmpSum; 
  
  for (int i = 0; i < this->NBodyValue; ++i)
    {
      for (int j = 0; j < this->NBodyValue - 1; ++j)
	{
	  int coef1 = (momFactor[j] / this->NbrLzValue) - (momFactor[j] / this->NbrLzValue) % this->NBodyValue;   
	  TmpIndices[j] = (double) (((i + momentumTransfer) % this->NBodyValue) - coef1);	
	  countIter[j] = 0;
	}
      Coefficient[i] = this->EvaluateGaussianSum(0, TmpIndices, 0, countIter, momFactor);
    }
  delete[] TmpIndices;
  delete[] momFactor;
  delete[] countIter;
  
  return Coefficient;
}



// evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
//
// nBodyValue = current index being incremented 
//TmpIndices = array containing the current value of all indices
// Sum = current value of the sum
// countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
// momFactor = array of indices that contains the information of the creation (or annihilation) indices
//return value = value of the coefficient
double ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor)
{
  if (nBodyValue == this->NBodyValue - 1)
    return 0.0;
  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double NormalizationCoefficient = pow(((double) this->NBodyValue) * M_PI * DoubleNbrLzValue,((double) (this->NBodyValue + 1)) / 4.0);
//  NormalizationCoefficient = 1.0;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Factor = 2.0*M_PI*this->Ratio / (DoubleNbrLzValue * ((double)(this->NBodyValue))* ((double)(this->NBodyValue)));
  int MinIter = 3;
  int TmpIndex = TmpIndices[0];
  double ExpFactor;

  
  for (int i = nBodyValue; i < this->NBodyValue - 1; ++i)
    {
      TmpIndices[i] = TmpIndex;
      countIter[i] = 0;
    }
  ExpFactor = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  double Coefficient = NormalizationCoefficient * exp(-Factor*ExpFactor);
  
  while ((Coefficient + Sum != Sum) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] += this->NBodyValue;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    ExpFactor = 0;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    Coefficient = NormalizationCoefficient * exp(-Factor*ExpFactor);
    
  }
  
  TmpIndices[nBodyValue] = TmpIndex - this->NBodyValue;
  countIter[nBodyValue] = 0;
  for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
  {
    TmpIndices[i] = TmpIndex;
    countIter[i] = 0;
  }
    
  ExpFactor = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  Coefficient = NormalizationCoefficient * exp(-Factor*ExpFactor);
  
  
  while ((Coefficient + Sum != Sum) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] -= this->NBodyValue;
    ExpFactor = 0;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    Coefficient = NormalizationCoefficient * exp(-Factor*ExpFactor);    
  }
  
  return Sum;
  
}
  
// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
//
// creationCoefficients = array that contains the creation coefficients
// annihilationCoefficients = array that contains the annihilation coefficients
// nbrPermutations1 = number of permutations of the creation indexes
// nbrPermutations2 = number of permutations of the annihilation indexes
// return value = numerical coefficient  

double ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients)
{
  double Coefficient = 0.0;
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);

 
  for (int j = 0; j < this->NBodyValue; ++j)
    Coefficient += creationCoefficients[j]*annihilationCoefficients[j];

  for (int i = 0; i < this->NBodyValue; ++i)
    Coefficient /= (M_PI * DoubleNbrLzValue);
  
  return (Coefficient/ (FactorialNBody.GetNumericalValue()));
}

// find the canonical index corresponding to a giving index
//
// index = linearized index 
// sign = reference on an optional sign resulting when finiding the canonical index
// indices = temporary array to store the full index
// nbrPermutations  = temporary array to store the numebr of permutations
// return value = canonical index

int ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::GetCanonicalIndex (int index, int& sign, int* indices, int* nbrPermutations)
{
  sign = 1;
  int g;
  int MomentumTransfer;
  this->GetIndicesFromLinearizedIndex (index, g, MomentumTransfer, indices);
  int TmpIndex = index;  
  int CountModulo;
  int TmpSign;
  nbrPermutations[0] = 0;
  SortArrayDownOrderingPermutationBubbleSort<int>(indices, this->NBodyValue, nbrPermutations[0]);
  int TmpSumGMomentumTransfer = (g + MomentumTransfer + this->NBodyValue) % this->NBodyValue;
  int SumGMomentumTransfer = TmpSumGMomentumTransfer;
  MomentumTransfer = - this->NBodyValue + 1;
  g = (TmpSumGMomentumTransfer + this->NBodyValue - 1) % this->NBodyValue;
  int CanonicalIndex = this->EvaluateLinearizedIndex(indices, g, MomentumTransfer);
  int IndexForCanonical = 0;
  for (int j = 0; j <  this->NbrLzValue - 1; ++j)
    {
      nbrPermutations[j + 1] = nbrPermutations[j];
      CountModulo = 0;
      for (int i = 0; i < this->NBodyValue; ++i)
	{
	  indices[i]++;
	  if (indices[i] == this->NbrLzValue)
	    {
	      CountModulo++;
	      indices[i] = 0;
	    }   
	}      
      
      if (CountModulo != 0)
	{
	  SortArrayDownOrderingPermutationBubbleSort<int>(indices, this->NBodyValue, nbrPermutations[j + 1]);	  
	}
      
      TmpSumGMomentumTransfer = (TmpSumGMomentumTransfer - CountModulo) % this->NBodyValue;
      MomentumTransfer = - this->NBodyValue + 1;
      g = (TmpSumGMomentumTransfer + this->NBodyValue - 1) % this->NBodyValue;
      TmpIndex = this->EvaluateLinearizedIndex(indices, g, MomentumTransfer);
      if (TmpIndex < CanonicalIndex)
	{
	  CanonicalIndex = TmpIndex;
	  IndexForCanonical = j + 1;
	}
    }
  
     
  if ((nbrPermutations[IndexForCanonical] % 2) == 0)
    sign = 1;
  else
    sign = -1;

  return CanonicalIndex;   
}


// test if creation/annihilation are in canonical form
//
// mIndices = array that contains the creation indices (will be modified)
// nIndices = array that contains the annihilation indices (will be modified)
// linearizedMIndices = linearized creation index
// linearizedNIndices = linearized annihilation index
// totalMomentum = momentum sector of the creation/annihilation indices
// canonicalMIndices = reference on the linearized creation index of the canonical form
// canonicalMIndices = reference on the linearized annihilation index of the canonical form
// canonicalTotalMomentum = reference on the momentum sector of the creation/annihilation indices of the canonical form
// canonicalPermutationCoefficient = additional permutation coefficient to get the canonical form
// return value = true if the indices are in canonical form

bool ParticleOnTorusWithSpinAndMagneticTranslationsNBodyHardCoreHamiltonian::FindCanonicalIndices(int* mIndices, int* nIndices, int linearizedMIndices, int linearizedNIndices, 
											   int totalMomentum, int& canonicalMIndices, int& canonicalNIndices, 
											   int& canonicalTotalMomentum, double& canonicalPermutationCoefficient)
{
  int* TmpMIndices = new int [this->NBodyValue];
  int* TmpNIndices = new int [this->NBodyValue];
  canonicalPermutationCoefficient = 1.0;
  canonicalMIndices = linearizedMIndices;
  canonicalNIndices = linearizedNIndices;
  canonicalTotalMomentum = totalMomentum;
  canonicalPermutationCoefficient = 1.0;
  bool IsCanonical = true;
  for (int k = 0; k < this->NBodyValue; ++k)
    {
      if (mIndices[k] > 0)
	{
	  TmpMIndices[k] = this->MaxMomentum - mIndices[k];
	}
      else
	{
	  TmpMIndices[k] = 0;
	}
      if (nIndices[k] > 0)
	{
	  TmpNIndices[k] = this->MaxMomentum - nIndices[k];
	}
      else
	{
	  TmpNIndices[k] = 0;
	}
    }
  int TmpNbrPermutation = 0;
  SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
  SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
  int TmpCanonicalTotalMomentum = 0;
  int TmpCanonicalMIndices = 0;
  int TmpCanonicalNIndices = 0;
  for (int k = this->NBodyValue - 1; k >= 0; --k)
    {
      TmpCanonicalMIndices *= this->MaxMomentum;
      TmpCanonicalMIndices += TmpMIndices[k];
      TmpCanonicalNIndices *= this->MaxMomentum;
      TmpCanonicalNIndices += TmpNIndices[k];
      TmpCanonicalTotalMomentum += TmpNIndices[k];
    }
  TmpCanonicalTotalMomentum %= this->MaxMomentum;
  if (TmpCanonicalMIndices > TmpCanonicalNIndices)
    {
      int Tmp = TmpCanonicalMIndices;
      TmpCanonicalMIndices = TmpCanonicalNIndices;
      TmpCanonicalNIndices = Tmp;
    }
  if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
      ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								 ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
    {
      IsCanonical = false;
      canonicalMIndices = TmpCanonicalMIndices;
      canonicalNIndices = TmpCanonicalNIndices;
      canonicalTotalMomentum = TmpCanonicalTotalMomentum;
      if ((TmpNbrPermutation & 1) == 0)
	{
	  canonicalPermutationCoefficient = 1.0;
	}
      else
	{
	  canonicalPermutationCoefficient = -1.0;
	}
    }
  for (int i = 1; i <  this->MaxMomentum; ++i)
    {
      bool ResuffleMFlag = false;
      bool ResuffleNFlag = false;
      for (int k = 0; k < this->NBodyValue; ++k)
	{
	  TmpMIndices[k] = mIndices[k] - i;
	  if (TmpMIndices[k] < 0)
	    {
	      TmpMIndices[k] += this->MaxMomentum;
	      ResuffleMFlag = true;
	    }
	  TmpNIndices[k] = nIndices[k] - i;
	  if (TmpNIndices[k] < 0)
	    {
	      TmpNIndices[k] += this->MaxMomentum;
	      ResuffleNFlag = true;
	    }
	}
      TmpNbrPermutation = 0;
      if (ResuffleMFlag == true)
	SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
      if (ResuffleNFlag == true)
	SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
      TmpCanonicalTotalMomentum = 0;
      TmpCanonicalMIndices = 0;
      TmpCanonicalNIndices = 0;
      for (int k = this->NBodyValue - 1; k >= 0; --k)
	{
	  TmpCanonicalMIndices *= this->MaxMomentum;
	  TmpCanonicalMIndices += TmpMIndices[k];
	  TmpCanonicalNIndices *= this->MaxMomentum;
	  TmpCanonicalNIndices += TmpNIndices[k];
	  TmpCanonicalTotalMomentum += TmpNIndices[k];
	}
      TmpCanonicalTotalMomentum %= this->MaxMomentum;
      if (TmpCanonicalMIndices > TmpCanonicalNIndices)
	{
	  int Tmp = TmpCanonicalMIndices;
	  TmpCanonicalMIndices = TmpCanonicalNIndices;
	  TmpCanonicalNIndices = Tmp;
	}
      if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
	  ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								     ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
	{
	  IsCanonical = false;
	  canonicalMIndices = TmpCanonicalMIndices;
	  canonicalNIndices = TmpCanonicalNIndices;
	  canonicalTotalMomentum = TmpCanonicalTotalMomentum;
	  if ((TmpNbrPermutation & 1) == 0)
	    {
	      canonicalPermutationCoefficient = 1.0;
	    }
	  else
	    {
	      canonicalPermutationCoefficient = -1.0;
	    }
 	}
      for (int k = 0; k < this->NBodyValue; ++k)
	{
	  if (TmpMIndices[k] > 0)
	    {
	      TmpMIndices[k] = this->MaxMomentum - TmpMIndices[k];
	    }
	  else
	    {
	      TmpMIndices[k] = 0;
	    }
	  if (TmpNIndices[k] > 0)
	    {
	      TmpNIndices[k] = this->MaxMomentum - TmpNIndices[k];
	    }
	  else
	    {
	      TmpNIndices[k] = 0;
	    }
	}
      // warning do not set TmpNbrPermutation to zero
      SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
      SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
      TmpCanonicalTotalMomentum = 0;
      TmpCanonicalMIndices = 0;
      TmpCanonicalNIndices = 0;
      for (int k = this->NBodyValue - 1; k >= 0; --k)
	{
	  TmpCanonicalMIndices *= this->MaxMomentum;
	  TmpCanonicalMIndices += TmpMIndices[k];
	  TmpCanonicalNIndices *= this->MaxMomentum;
	  TmpCanonicalNIndices += TmpNIndices[k];
	  TmpCanonicalTotalMomentum += TmpNIndices[k];
	}
      TmpCanonicalTotalMomentum %= this->MaxMomentum;
      if (TmpCanonicalMIndices > TmpCanonicalNIndices)
	{
	  int Tmp = TmpCanonicalMIndices;
	  TmpCanonicalMIndices = TmpCanonicalNIndices;
	  TmpCanonicalNIndices = Tmp;
	}
      if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
	  ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								     ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
	{
	  IsCanonical = false;
	  canonicalMIndices = TmpCanonicalMIndices;
	  canonicalNIndices = TmpCanonicalNIndices;
	  canonicalTotalMomentum = TmpCanonicalTotalMomentum;
	  if ((TmpNbrPermutation & 1) == 0)
	    {
	      canonicalPermutationCoefficient = 1.0;
	    }
	  else
	    {
	      canonicalPermutationCoefficient = -1.0;
	    }
	}
    }
  delete[] TmpMIndices;
  delete[] TmpNIndices;
  return IsCanonical;
}

