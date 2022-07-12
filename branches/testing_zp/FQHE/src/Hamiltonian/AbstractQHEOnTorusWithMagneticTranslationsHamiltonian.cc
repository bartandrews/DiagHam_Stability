////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                                                                            //
//                        last modification : 18/11/2003                      //
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
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::hex;
using std::dec;
using std::endl;
using std::ostream;


// destructor
//

AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::~AbstractQHEOnTorusWithMagneticTranslationsHamiltonian()
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Cosinus;
  double Sinus;
  int NbrTranslation;
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  Complex TmpZ;
  Complex Z (0.0, 0.0);
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
    {
      m1 = this->M1Value[j];
      m2 = this->M2Value[j];
      m3 = this->M3Value[j];
      m4 = this->M4Value[j];
      TmpInteraction = this->InteractionFactors[j];
      for (int i = 0; i < Dim; ++i)
	{
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < Dim)
	    {
	      Coefficient *= TmpInteraction;
	      Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
	      Sinus = Coefficient * this->SinusTable[NbrTranslation];
	      TmpZ.Re = ((V2.Re(i) * Cosinus) - (V2.Im(i) * Sinus));
	      TmpZ.Im += ((V2.Re(i) * Sinus) + (V2.Im(i) * Cosinus));
	      Z += Conj(V1[Index]) * TmpZ;
	    }
	}
    }
  m1 = this->M1Value[ReducedNbrInteractionFactors];
  m2 = this->M2Value[ReducedNbrInteractionFactors];
  m3 = this->M3Value[ReducedNbrInteractionFactors];
  m4 = this->M4Value[ReducedNbrInteractionFactors];
  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
      if (Index < Dim)
	{
	  Coefficient *= TmpInteraction;
	  Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
	  Sinus = Coefficient * this->SinusTable[NbrTranslation];
	  TmpZ.Re = ((V2.Re(i) * Cosinus) - (V2.Im(i) * Sinus));
	  TmpZ.Im += ((V2.Re(i) * Sinus) + (V2.Im(i) * Cosinus));
	  Z += Conj(V1[Index]) * TmpZ;
	}
      Z += this->EnergyShift * (Conj(V1[i]) * V2[i]);
    }
  return Complex(Z);
}
  
// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
											       int firstComponent, int nbrComponent) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
												  int firstComponent, int nbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
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

ComplexVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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
ComplexVector& AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = 0.0;
  double Coefficient;
  double Cosinus;
  double Sinus;
  int NbrTranslation;
  if (this->FastMultiplicationFlag == false)
    {
      double Coefficient;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	      if (Index < Dim)
		{
		  Coefficient *= TmpInteraction;
		  Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
		  Sinus = Coefficient * this->SinusTable[NbrTranslation];
		  vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
		  vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
		}
	    }
	}
      m1 = this->M1Value[ReducedNbrInteractionFactors];
      m2 = this->M2Value[ReducedNbrInteractionFactors];
      m3 = this->M3Value[ReducedNbrInteractionFactors];
      m4 = this->M4Value[ReducedNbrInteractionFactors];
      TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < Dim)
	    {
	      Coefficient *= TmpInteraction;
	      Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
	      Sinus = Coefficient * this->SinusTable[NbrTranslation];
	      vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
	      vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
	    }
	  vDestination.Re(i) += this->EnergyShift * vSource.Re(i);
	  vDestination.Im(i) += this->EnergyShift * vSource.Im(i); 
	}
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int* TmpNbrTranslationArray;
	  int j;
	  int TmpNbrInteraction;
	  double TmpRe;
	  double TmpIm;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[i];
	      TmpRe = vSource.Re(i);
	      TmpIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Cosinus = TmpCoefficientArray[j];
		  NbrTranslation = TmpNbrTranslationArray[j];
		  Sinus = Cosinus * this->SinusTable[NbrTranslation];
		  Cosinus *= this->CosinusTable[NbrTranslation];
		  vDestination.Re(TmpIndexArray[j]) += ((Cosinus * TmpRe) - (Sinus * TmpIm));
		  vDestination.Im(TmpIndexArray[j]) += ((Sinus * TmpRe) + (Cosinus * TmpIm));
		}
	      vDestination.Re(i) += Shift * TmpRe;
	      vDestination.Im(i) += Shift * TmpIm;
	    }
	}
      else
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int* TmpNbrTranslationArray;
	  int j;
	  int TmpNbrInteraction;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  double TmpRe;
	  double TmpIm;
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
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[Pos];
	      TmpRe = vSource.Re(i);
	      TmpIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Cosinus = TmpCoefficientArray[j];
		  NbrTranslation = TmpNbrTranslationArray[j];
		  Sinus = Cosinus * this->SinusTable[NbrTranslation];
		  Cosinus *= this->CosinusTable[NbrTranslation];
		  vDestination.Re(TmpIndexArray[j]) += ((Cosinus * TmpRe) - (Sinus * TmpIm));
		  vDestination.Im(TmpIndexArray[j]) += ((Sinus * TmpRe) + (Cosinus * TmpIm));
		}
	      vDestination.Re(i) += Shift * TmpRe;
	      vDestination.Im(i) += Shift * TmpIm;
	      ++Pos;
	    }
	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
	  double TmpInteraction;
	  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
		  {
		    m1 = this->M1Value[j];
		    m2 = this->M2Value[j];
		    m3 = this->M3Value[j];
		    m4 = this->M4Value[j];
		    TmpInteraction = this->InteractionFactors[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
			Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
			if (Index < Dim)
			  {
			    Coefficient *= TmpInteraction;
			    Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
			    Sinus = Coefficient * this->SinusTable[NbrTranslation];
			    vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
			    vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
			  }
		      }
		  }
		m1 = this->M1Value[ReducedNbrInteractionFactors];
		m2 = this->M2Value[ReducedNbrInteractionFactors];
		m3 = this->M3Value[ReducedNbrInteractionFactors];
		m4 = this->M4Value[ReducedNbrInteractionFactors];
		TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
		    if (Index < Dim)
		      {
			Coefficient *= TmpInteraction;
			Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
			Sinus = Coefficient * this->SinusTable[NbrTranslation];
			vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
			vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
		      }
		    vDestination.Re(i) += this->EnergyShift * vSource.Re(i);
		    vDestination.Im(i) += this->EnergyShift * vSource.Im(i);
		  }
	      }
	}
   }
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::FastMultiplicationMemory(long allowedMemory)
{

  this->NbrInteractionPerComponent = new int [this->Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension();
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension()) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
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

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  int NbrTranslation;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  int LastComponent = lastComponent + firstComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i];
	    }
	}    
    }
  Memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * this->Particles->GetHilbertSpaceDimension() + 
	    Memory *  (sizeof (int) + sizeof(double) + sizeof(int)));
  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::EnableFastMultiplication()
{
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  int NbrTranslation;
  double* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];
  this->InteractionPerComponentNbrTranslation = new int* [ReducedSpaceDimension];

  int TotalPos = 0;
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
      this->InteractionPerComponentNbrTranslation[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[TotalPos];
      Pos = 0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      TmpNbrTranslationArray[Pos] = NbrTranslation;
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

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  int NbrTranslation;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  int Pos;
  int Min = firstComponent / this->FastMultiplicationStep;
  int Max = lastComponent / this->FastMultiplicationStep;
  
  for (int i = Min; i < Max; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];      
      this->InteractionPerComponentNbrTranslation[i] = new int [this->NbrInteractionPerComponent[i]];
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[i];
      Pos = 0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  Index = this->Particles->AdAdAA(i * this->FastMultiplicationStep, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      TmpNbrTranslationArray[Pos] = NbrTranslation;
	      ++Pos;
	    }
	}
    }
}


// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::SavePrecalculation (char* fileName)
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
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentNbrTranslation[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
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

bool AbstractQHEOnTorusWithMagneticTranslationsHamiltonian::LoadPrecalculation (char* fileName)
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
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentNbrTranslation[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

