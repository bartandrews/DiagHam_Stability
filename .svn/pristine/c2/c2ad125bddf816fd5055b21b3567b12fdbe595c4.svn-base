////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          hardcore 3-body interaction                       //
//                                                                            //
//                        last modification : 30/04/2014                      //
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


#include "config.h"
#include "Hamiltonian/ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <iostream>


// default constructor
//

ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// twoBodyDeltaStrength = strength of the additional two body delta interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
																	 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
																	 char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->NBodyValue = 3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TwoBodyFlag = false;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  this->PrecalculatedInteractionCoefficients = new double***** [this->NbrLzValue];
  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
    {
      this->PrecalculatedInteractionCoefficients[m1] = new double**** [this->NbrLzValue];
      for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
	{
	  this->PrecalculatedInteractionCoefficients[m1][m2] = new double*** [this->NbrLzValue];
	  for (int m3 = 0; m3 < this->NbrLzValue; ++m3)
	    {
	      this->PrecalculatedInteractionCoefficients[m1][m2][m3] = new double** [this->NbrLzValue];
	      for (int n1 = 0; n1 < this->NbrLzValue; ++n1)
		{
		  this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1] = new double* [this->NbrLzValue];
		  for (int n2 = 0; n2 < this->NbrLzValue; ++n2)
		    {
		      this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2] = new double [this->NbrLzValue];
		    }
		}
	    }
	}
    }
  char* InteractionCoefficientFileName = new char [512];
  sprintf (InteractionCoefficientFileName, "threebodydelta_interactioncoefficient_2s_%d_ratio_%.10f.dat", this->NbrLzValue, this->Ratio);
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
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    {
	      for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
		{
		  for (int m3 = 0; m3 < this->NbrLzValue; ++m3)
		    {
		      for (int n1 = 0; n1 < this->NbrLzValue; ++n1)
			{
			  for (int n2 = 0; n2 < this->NbrLzValue; ++n2)
			    {
			      ReadBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2], this->NbrLzValue);
// 			      for (int n3 = 0; n3 < this->NbrLzValue; ++n3)
// 				{
// 				  ReadLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3]);
// 				}
			    }
			}
		    }
		}
	    }
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
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    {
	      for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
		{
		  for (int m3 = 0; m3 < this->NbrLzValue; ++m3)
		    {
		      for (int n1 = 0; n1 < this->NbrLzValue; ++n1)
			{
			  for (int n2 = 0; n2 < this->NbrLzValue; ++n2)
			    {
			      for (int n3 = 0; n3 < this->NbrLzValue; ++n3)
				{
				  this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3] = this->EvaluateInteractionCoefficient(m1, m2, m3, n1, n2, n3);
				}
			      WriteBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2], this->NbrLzValue);
			    }
			}
		    }
		}
	    }
	  File.close();
	}
    }
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
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

ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::~ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian()
{
}
  
// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// evaluate all interaction factors
//   

void ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0l;
  this->GetIndices();
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
    }
  else
    {
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      int n1 = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      int n2  = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      int n3  = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  int m1 = this->NBodySectorIndicesPerSum[i][j2 * 3];
		  int m2  = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 1];
		  int m3  = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 2];		  
		  double TmpInteraction = 0.0;
		  if (m1 != m2)
		    {
		      if (m2 != m3)
			{
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m2, m3, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m3, m2, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m2, m1, m3, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m2, m3, m1, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m3, m1, m2, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m3, m2, m1, n1, n2, n3);
			}
		      else
			{
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m2, m3, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m2, m1, m3, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m2, m3, m1, n1, n2, n3);
			}
		    }
		  else
		    {		      
		      if (m2 != m3)
			{
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m2, m3, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m3, m2, n1, n2, n3);
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m3, m1, m2, n1, n2, n3);
			}
		      else
			{
			  TmpInteraction += this->EvaluateInteractionNIndexSymmetrizedCoefficient(m1, m2, m3, n1, n2, n3);
			}
		    }
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction;
		  
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}
	
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term
//
// m1 = first creation operator index
// m2 = second creation operator index
// m3 = third creation operator index
// n1 = first annihilation operator index
// n2 = second annihilation operator index
// n3 = thrid annihilation operator index
//
// return value = numerical coefficient  

double ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int n1, int n2, int n3)
{
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double ResultNy2 = 0.0;
  double ResultNy1 = 0.0;
  double ResultNx2 = 0.0;
  double ResultNx1 = 0.0;
  double Nx1;
  double Nx2;
  double Q1;
  double Q2;
  double Q3;
  double Q4;
  double Q5;
  double IniNy1 = (double) (m1 - n3);
  double Ny1 = IniNy1;
  double IniNy2 = (double) (n1 - m3);
  double Ny2 = IniNy2;
  double PremFactor1 = ((2.0 * ((double) (m1 - n2))) - Ny2)* PIOnM;
  double PremFactor2 = ((2.0 * ((double) (n2 - m3))) - Ny1) * PIOnM;
  double Factor1 = PremFactor1;
  double Factor2 = PremFactor2;
  double Precision = 1.0;
  double Precision1 = 1.0;
  double Coefficient = 1.0;
  double Coefficient1 = 1.0;
  double Coefficient2 = 1.0;
  //cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << n1 << " "<< n2 << " "<< n3 << endl;
  
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      ResultNy1 = 0.0;
      Ny1 = IniNy1;
      Factor2 = PremFactor2;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  
	  
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2 = 0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 =  2.0 * cos(Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2 * ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 -= M_PI;
	  Ny1 += DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 -= M_PI;
      Ny2 += DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2 - DoubleNbrLzValue;
  Factor1 = PremFactor1 + M_PI;
  Factor2 = PremFactor2;
  Coefficient = 1.0;
  
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      ResultNy1 = 0.0;
      Ny1 = IniNy1;
      Factor2 = PremFactor2;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2* ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1;
	  Factor2 -= M_PI;
	  Ny1 += DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 += M_PI;
      Ny2 -= DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2;
  Factor1 = PremFactor1;
  Factor2 = PremFactor2 + M_PI;
  
  Coefficient = 1.0;
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      
      Ny1 = IniNy1 - DoubleNbrLzValue;
      Factor2 = PremFactor2 + M_PI;
      ResultNy1 = 0.0;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2 * ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 += M_PI;
	  Ny1 -= DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 -= M_PI;
      Ny2 += DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2 - DoubleNbrLzValue;
  Factor1 = PremFactor1 + M_PI;
  Factor2 = PremFactor2 + M_PI;
  
  Coefficient = 1.0;	
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      
      Ny1 = IniNy1 - DoubleNbrLzValue;
      Factor2 = PremFactor2 + M_PI;
      ResultNy1 = 0.0;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
			}
	      ResultNx1 += Coefficient2* ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 += M_PI;
	  Ny1 -= DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 += M_PI;
      Ny2 -= DoubleNbrLzValue;
    }
  
  return (ResultNy2 / (24.0 * (M_PI * DoubleNbrLzValue)*(M_PI * DoubleNbrLzValue)));
}
  
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  return 0.0;
}

