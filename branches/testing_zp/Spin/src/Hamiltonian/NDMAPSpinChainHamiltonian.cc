////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin chain hamiltonian with translations             //
//                                                                            //
//                        last modification : 04/03/2002                      //
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


#include "Hamiltonian/NDMAPSpinChainHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/NDMAPPrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// chain = pointer to Hilbert space of the associated system
// nbrSpin = number of spin
// j = coupling constants between spins
// jz = coupling constants between spins in the z direction
// parallelMagneticField = magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction parallel to the chain
// perpendicularMagneticField = magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction perpendicular to the chain
// d = single ion anisotropy constant
// e = single ion in-plane anisotropy constant
// architecture = architecture to use for precalculation
// memory = amount of memory (in bytes) that can allocated for precalcutation

NDMAPSpinChainHamiltonian::NDMAPSpinChainHamiltonian(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double jz, double parallelMagneticField, 
						     double perpendicularMagneticField, double d, double e,  AbstractArchitecture* architecture, unsigned long memory)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = j;
  this->HalfJ = this->J * 0.5;
  this->Jz = jz;
  this->ParallelMagneticField = parallelMagneticField;
  this->PerpendicularMagneticField = perpendicularMagneticField;
  this->HalfPerpendicularMagneticField = 0.5 * this->PerpendicularMagneticField;
  this->D = d;
  this->E = e;
  this->HalfE = 0.5 * this->E;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;
  if ((memory > 0) && (this->Architecture != 0))
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

// destructor
//

NDMAPSpinChainHamiltonian::~NDMAPSpinChainHamiltonian() 
{
  delete[] this->SzSzContributions;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  if (  this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Chain->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	  delete[] this->InteractionPerComponentNbrTranslations[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      delete[] this->InteractionPerComponentNbrTranslations;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void NDMAPSpinChainHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* NDMAPSpinChainHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

NDMAPSpinChainHamiltonian& NDMAPSpinChainHamiltonian::SetChain(AbstractSpinChainWithTranslations* chain)
{  
  delete[] this->SzSzContributions;
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  return *this;  
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int NDMAPSpinChainHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void NDMAPSpinChainHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex NDMAPSpinChainHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex(0);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex NDMAPSpinChainHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  Complex Z (0.0, 0.0);
  Complex TmpZ;
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
    {
      Z += this->SzSzContributions[i] * (Conj(V1[i]) * V2[i]);
    }
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				(V2.Im(i) * this->SinusTable[NbrTranslation]));
	      TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				(V2.Im(i) * this->CosinusTable[NbrTranslation]));
	      Z += Conj(V1[pos]) * TmpZ;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
				(V2.Im(i) * this->SinusTable[NbrTranslation]));
	      TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
				(V2.Im(i) * this->CosinusTable[NbrTranslation]));
	      Z += Conj(V1[pos]) * TmpZ;
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  TmpZ.Re =Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
			   (V2.Im(i) * this->SinusTable[NbrTranslation]));
	  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
			    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
	  Z += Conj(V1[pos]) * TmpZ;
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
			    (V2.Im(i) * this->SinusTable[NbrTranslation]));
	  TmpZ.Im = Coef * ((V2.Re(i) * this->SinusTable[NbrTranslation]) +
			    (V2.Im(i) * this->CosinusTable[NbrTranslation]));
	  Z += Conj(V1[pos]) * TmpZ;
	}      
    }
  if (this->PerpendicularMagneticField != 0.0)
    {
      for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfPerpendicularMagneticField;
		  TmpZ.Re = Coef * ((this->CosinusTable[NbrTranslation] * V2.Re(i))
				    - (this->SinusTable[NbrTranslation] * V2.Im(i)));
		  TmpZ.Im = Coef * ((this->SinusTable[NbrTranslation] * V2.Re(i))
				    + (this->CosinusTable[NbrTranslation] * V2.Im(i)));
		  Z += Conj(V1[pos]) * TmpZ;
		}
	      pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfPerpendicularMagneticField;
		  TmpZ.Re = Coef * ((this->CosinusTable[NbrTranslation] * V2.Re(i))
				    - (this->SinusTable[NbrTranslation] * V2.Im(i)));
		  TmpZ.Im = Coef * ((this->SinusTable[NbrTranslation] * V2.Re(i))
				    + (this->CosinusTable[NbrTranslation] * V2.Im(i)));
		  Z += Conj(V1[pos]) * TmpZ;
		}
	    }
	}
    }
  if (this->E != 0.0)
    {
      for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfE;
		  TmpZ.Re = Coef * ((this->CosinusTable[NbrTranslation] * V2.Re(i))
				    - (this->SinusTable[NbrTranslation] * V2.Im(i)));
		  TmpZ.Im = Coef * ((this->SinusTable[NbrTranslation] * V2.Re(i))
				    + (this->CosinusTable[NbrTranslation] * V2.Im(i)));
		  Z += Conj(V1[pos]) * TmpZ;
		}
	      pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfE;
		  TmpZ.Re = Coef * ((this->CosinusTable[NbrTranslation] * V2.Re(i))
				    - (this->SinusTable[NbrTranslation] * V2.Im(i)));
		  TmpZ.Im = Coef * ((this->SinusTable[NbrTranslation] * V2.Re(i))
				    + (this->CosinusTable[NbrTranslation] * V2.Im(i)));
		  Z += Conj(V1[pos]) * TmpZ;
		}      
	    }
	}
    }
  return Z;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& NDMAPSpinChainHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      for (int i = firstComponent; i < Last; i++)
	{
	  vDestination.Re(i) = this->SzSzContributions[i] * vSource.Re(i);
	  vDestination.Im(i) = this->SzSzContributions[i] * vSource.Im(i);
	}
      for (int i = firstComponent; i < Last; i++)
	{
	  for (int j = 0; j < MaxPos; j++)
	    {
	      pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfJ;
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	      pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfJ;
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	    }    
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	  pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }      
	}
      if (this->PerpendicularMagneticField != 0.0)
	{
	  for (int i = firstComponent; i < Last; i++)
	    {
	      for (int j = 0; j < this->NbrSpin; j++)
		{
		  pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfPerpendicularMagneticField;
		      vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
						      - (this->SinusTable[NbrTranslation] * vSource.Im(i)));
		      vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
						      + (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
		    }
		  pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfPerpendicularMagneticField;
		      vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
						      - (this->SinusTable[NbrTranslation] * vSource.Im(i)));
		      vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
						      + (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
		    }
		}
	    }
	}
      if (this->E != 0.0)
	{
	  for (int i = firstComponent; i < Last; i++)
	    {
	      for (int j = 0; j < this->NbrSpin; j++)
		{
		  pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfE;
		      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		    }
		  pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfE;
		      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		    }      
		}
	    }
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
	  double CoefficientRe;
	  double CoefficientIm;
	  for (int i = firstComponent; i < Last; i++)
	    {
	      vDestination.Re(i) = this->SzSzContributions[i] * vSource.Re(i);
	      vDestination.Im(i) = this->SzSzContributions[i] * vSource.Im(i);
	    }
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslations[i];
	      CoefficientRe = vSource.Re(i);
	      CoefficientIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination.Re(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->CosinusTable[TmpNbrTranslationArray[j]] - 
										  CoefficientIm * this->SinusTable[TmpNbrTranslationArray[j]]);
		  vDestination.Im(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->SinusTable[TmpNbrTranslationArray[j]] + 
										  CoefficientIm * this->CosinusTable[TmpNbrTranslationArray[j]]);
		}
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
	  double CoefficientRe;
	  double CoefficientIm;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = firstComponent; i < Last; i++)
	    {
	      vDestination.Re(i) = this->SzSzContributions[i] * vSource.Re(i);
	      vDestination.Im(i) = this->SzSzContributions[i] * vSource.Im(i);
	    }
	  for (int i = PosMod + firstComponent; i < Last; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslations[Pos];
	      CoefficientRe = vSource.Re(i);
	      CoefficientIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination.Re(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->CosinusTable[TmpNbrTranslationArray[j]] - 
										  CoefficientIm * this->SinusTable[TmpNbrTranslationArray[j]]);
		  vDestination.Im(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->SinusTable[TmpNbrTranslationArray[j]] + 
										  CoefficientIm * this->CosinusTable[TmpNbrTranslationArray[j]]);
		}
	      ++Pos;
	    }
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {	
		for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		  {
		    for (int j = 0; j < MaxPos; j++)
		      {
			pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
			if (pos != this->Chain->GetHilbertSpaceDimension())
			  {
			    Coef *= this->HalfJ;
			    vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							    (vSource.Im(i) * this->SinusTable[NbrTranslation]));
			    vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							    (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			  }
			pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
			if (pos != this->Chain->GetHilbertSpaceDimension())
			  {
			    Coef *= this->HalfJ;
			    vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							    (vSource.Im(i) * this->SinusTable[NbrTranslation]));
			    vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							    (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			  }
		      }    
		    pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
		    if (pos != this->Chain->GetHilbertSpaceDimension())
		      {
			Coef *= this->HalfJ;
			vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							(vSource.Im(i) * this->SinusTable[NbrTranslation]));
			vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		      }
		    pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
		    if (pos != this->Chain->GetHilbertSpaceDimension())
		      {
			Coef *= this->HalfJ;
			vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							(vSource.Im(i) * this->SinusTable[NbrTranslation]));
			vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		      }      
		  }
		if (this->PerpendicularMagneticField != 0.0)
		  {
		    for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		      {
			for (int j = 0; j < this->NbrSpin; j++)
			  {
			    pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfPerpendicularMagneticField;
				vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
								- (this->SinusTable[NbrTranslation] * vSource.Im(i)));
				vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
								+ (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
			      }
			    pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfPerpendicularMagneticField;
				vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
								- (this->SinusTable[NbrTranslation] * vSource.Im(i)));
				vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
								+ (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
			      }
			  }
		      }
		  }
		if (this->E != 0.0)
		  {
		    for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		      {
			for (int j = 0; j < this->NbrSpin; j++)
			  {
			    pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfE;
				vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
								(vSource.Im(i) * this->SinusTable[NbrTranslation]));
				vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
								(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			      }
			    pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfE;
				vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
								(vSource.Im(i) * this->SinusTable[NbrTranslation]));
				vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
								(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			      }      
			  }
		      }
		  }
	      }	
	}
    }
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

ComplexVector& NDMAPSpinChainHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      for (int i = firstComponent; i < Last; i++)
	{
	  vDestination.Re(i) += this->SzSzContributions[i] * vSource.Re(i);
	  vDestination.Im(i) += this->SzSzContributions[i] * vSource.Im(i);
	  for (int j = 0; j < MaxPos; j++)
	    {
	      pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfJ;
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	      pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  Coef *= this->HalfJ;
		  vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		  vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		}
	    }    
	  pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	  pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
					      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
	      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
					      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
	    }      
	}
      if (this->PerpendicularMagneticField != 0.0)
	{
	  for (int i = firstComponent; i < Last; i++)
	    {
	      for (int j = 0; j < this->NbrSpin; j++)
		{
		  pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfPerpendicularMagneticField;
		      vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
						      - (this->SinusTable[NbrTranslation] * vSource.Im(i)));
		      vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
						      + (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
		    }
		  pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfPerpendicularMagneticField;
		      vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
						      - (this->SinusTable[NbrTranslation] * vSource.Im(i)));
		      vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
						      + (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
		    }
		}
	    }
	}
      if (this->E != 0.0)
	{
	  for (int i = firstComponent; i < Last; i++)
	    {
	      for (int j = 0; j < this->NbrSpin; j++)
		{
		  pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfE;
		      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		    }
		  pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
		  if (pos != this->Chain->GetHilbertSpaceDimension())
		    {
		      Coef *= this->HalfE;
		      vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSource.Im(i) * this->SinusTable[NbrTranslation]));
		      vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		    }      
		}
	    }
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
	  double CoefficientRe;
	  double CoefficientIm;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslations[i];
	      CoefficientRe = vSource.Re(i);
	      CoefficientIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination.Re(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->CosinusTable[TmpNbrTranslationArray[j]] - 
										  CoefficientIm * this->SinusTable[TmpNbrTranslationArray[j]]);
		  vDestination.Im(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->SinusTable[TmpNbrTranslationArray[j]] + 
										  CoefficientIm * this->CosinusTable[TmpNbrTranslationArray[j]]);
		}
	      vDestination.Re(i) +=  this->SzSzContributions[i] * CoefficientRe;
	      vDestination.Im(i) +=  this->SzSzContributions[i] * CoefficientIm;
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
	  double CoefficientRe;
	  double CoefficientIm;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = PosMod + firstComponent; i < Last; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslations[Pos];
	      CoefficientRe = vSource.Re(i);
	      CoefficientIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination.Re(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->CosinusTable[TmpNbrTranslationArray[j]] - 
										  CoefficientIm * this->SinusTable[TmpNbrTranslationArray[j]]);
		  vDestination.Im(TmpIndexArray[j]) +=  TmpCoefficientArray[j] * (CoefficientRe * this->SinusTable[TmpNbrTranslationArray[j]] + 
										  CoefficientIm * this->CosinusTable[TmpNbrTranslationArray[j]]);
		}
	      vDestination.Re(i) +=  this->SzSzContributions[i] * CoefficientRe;
	      vDestination.Im(i) +=  this->SzSzContributions[i] * CoefficientIm;
	      ++Pos;
	    }
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {	
		for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		  {
		    vDestination.Re(i) += this->SzSzContributions[i] * vSource.Re(i);
		    vDestination.Im(i) += this->SzSzContributions[i] * vSource.Im(i);
		    for (int j = 0; j < MaxPos; j++)
		      {
			pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
			if (pos != this->Chain->GetHilbertSpaceDimension())
			  {
			    Coef *= this->HalfJ;
			    vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							    (vSource.Im(i) * this->SinusTable[NbrTranslation]));
			    vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							    (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			  }
			pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
			if (pos != this->Chain->GetHilbertSpaceDimension())
			  {
			    Coef *= this->HalfJ;
			    vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							    (vSource.Im(i) * this->SinusTable[NbrTranslation]));
			    vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							    (vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			  }
		      }    
		    pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
		    if (pos != this->Chain->GetHilbertSpaceDimension())
		      {
			Coef *= this->HalfJ;
			vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							(vSource.Im(i) * this->SinusTable[NbrTranslation]));
			vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		      }
		    pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
		    if (pos != this->Chain->GetHilbertSpaceDimension())
		      {
			Coef *= this->HalfJ;
			vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
							(vSource.Im(i) * this->SinusTable[NbrTranslation]));
			vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
							(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
		      }      
		  }
		if (this->PerpendicularMagneticField != 0.0)
		  {
		    for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		      {
			for (int j = 0; j < this->NbrSpin; j++)
			  {
			    pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfPerpendicularMagneticField;
				vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
								- (this->SinusTable[NbrTranslation] * vSource.Im(i)));
				vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
								+ (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
			      }
			    pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfPerpendicularMagneticField;
				vDestination.Re(pos) += Coef * ((this->CosinusTable[NbrTranslation] * vSource.Re(i))
								- (this->SinusTable[NbrTranslation] * vSource.Im(i)));
				vDestination.Im(pos) += Coef * ((this->SinusTable[NbrTranslation] * vSource.Re(i))
								+ (this->CosinusTable[NbrTranslation] * vSource.Im(i)));
			      }
			  }
		      }
		  }
		if (this->E != 0.0)
		  {
		    for (int i = firstComponent + k; i < Last; i += this->FastMultiplicationStep)
		      {
			for (int j = 0; j < this->NbrSpin; j++)
			  {
			    pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfE;
				vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
								(vSource.Im(i) * this->SinusTable[NbrTranslation]));
				vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
								(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			      }
			    pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
			    if (pos != this->Chain->GetHilbertSpaceDimension())
			      {
				Coef *= this->HalfE;
				vDestination.Re(pos) += Coef * ((vSource.Re(i) * this->CosinusTable[NbrTranslation]) -
								(vSource.Im(i) * this->SinusTable[NbrTranslation]));
				vDestination.Im(pos) += Coef * ((vSource.Re(i) * this->SinusTable[NbrTranslation]) +
								(vSource.Im(i) * this->CosinusTable[NbrTranslation]));
			      }      
			  }
		      }
		  }
	      }	
	}
    }
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> NDMAPSpinChainHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> NDMAPSpinChainHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void NDMAPSpinChainHamiltonian::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
    }
}

// evaluate diagonal matrix elements
// 

void NDMAPSpinChainHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double HalfMagneticField = 0.5 * this->ParallelMagneticField;
  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
      this->SzSzContributions[i] *= this->Jz;
      this->SzSzContributions[i] += HalfMagneticField * ((double) this->Chain->TotalSz(i));
      this->SzSzContributions[i] += this->D * this->Chain->TotalSzSz(i);      
    }
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long NDMAPSpinChainHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  this->NbrInteractionPerComponent = new int [this->Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
 
  NDMAPPrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - ((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * this->Chain->GetHilbertSpaceDimension();
  if ((TmpMemory < 0) || ((TmpMemory / ((int) ((2 * sizeof (int)) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) ((2 * sizeof (int)) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
	  if (this->Chain->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - ((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
      Memory = (((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
    }
  else
    {
      Memory = ((((2 * sizeof (int*)) + sizeof (int) + sizeof(double*)) * this->Chain->GetHilbertSpaceDimension()) + 
		(Memory * ((2 * sizeof (int)) + sizeof(double))));
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

long NDMAPSpinChainHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int pos;
  double Coef;
  long Memory = 0;
  int LastComponent = lastComponent + firstComponent;
  int MaxPos = this->NbrSpin - 1;
  int NbrTranslation;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < MaxPos; ++j) 
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i];
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i];
	    }
	}    
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  ++Memory;
	  ++this->NbrInteractionPerComponent[i];	  
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  ++Memory;
	  ++this->NbrInteractionPerComponent[i];
	}
    }
  if (this->PerpendicularMagneticField != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i];
		}
	      pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i];
		}
	    }
	}
    }
  if (this->E != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; i++)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i];
		}
	      pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i];
		}
	    }
	}
    }
  return Memory;
}

// enable fast multiplication algorithm
//

void NDMAPSpinChainHamiltonian::EnableFastMultiplication()
{
  int Index;
  int MaxPos = this->NbrSpin - 1;
  double Coef;
  int NbrTranslation;
  int* TmpIndexArray;
  int* TmpNbrTranslationArray;
  double* TmpCoefficientArray;
  int pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = this->Chain->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != this->Chain->GetHilbertSpaceDimension())
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];
  this->InteractionPerComponentNbrTranslations = new int* [ReducedSpaceDimension];

  int TotalPos = 0;
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
      this->InteractionPerComponentNbrTranslations[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslations[TotalPos];
      Index = 0;
      for (int j = 0; j < MaxPos; ++j) 
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Index] = pos;
	      TmpCoefficientArray[Index] = Coef * this->HalfJ;
	      TmpNbrTranslationArray[Index] = NbrTranslation;
	      ++Index;
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Index] = pos;
	      TmpCoefficientArray[Index] = Coef * this->HalfJ;
	      TmpNbrTranslationArray[Index] = NbrTranslation;
	      ++Index;
	    }
	}    
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  TmpIndexArray[Index] = pos;
	  TmpCoefficientArray[Index] = Coef * this->HalfJ;
	  TmpNbrTranslationArray[Index] = NbrTranslation;
	  ++Index;
	}
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  TmpIndexArray[Index] = pos;
	  TmpCoefficientArray[Index] = Coef * this->HalfJ;
	  TmpNbrTranslationArray[Index] = NbrTranslation;
	  ++Index;
	}
      if (this->PerpendicularMagneticField != 0.0)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->Spi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Index] = pos;
		  TmpCoefficientArray[Index] = Coef * this->HalfPerpendicularMagneticField;
		  TmpNbrTranslationArray[Index] = NbrTranslation;
		  ++Index;
		}
	      pos = this->Chain->Smi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Index] = pos;
		  TmpCoefficientArray[Index] = Coef * this->HalfPerpendicularMagneticField;
		  TmpNbrTranslationArray[Index] = NbrTranslation;
		  ++Index;
		}
	    }
	}
      if (this->E != 0.0)
	{
	  for (int j = 0; j < this->NbrSpin; j++)
	    {
	      pos = this->Chain->SpiSpi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Index] = pos;
		  TmpCoefficientArray[Index] = Coef * this->HalfE;
		  TmpNbrTranslationArray[Index] = NbrTranslation;
		  ++Index;
		}
	      pos = this->Chain->SmiSmi(j, i, Coef, NbrTranslation);
	      if (pos != this->Chain->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Index] = pos;
		  TmpCoefficientArray[Index] = Coef * this->HalfE;
		  TmpNbrTranslationArray[Index] = NbrTranslation;
		  ++Index;
		}
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

void NDMAPSpinChainHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, NDMAPSpinChainHamiltonian& H) 
{
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	  TmpV2.Im(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      TmpV2.Im(i) = 0.0;
      ((AbstractHamiltonian*) &H)->LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain->GetHilbertSpaceDimension(); j++)
	{
	  Str << "(" << TmpV[j].Re(i) << ", " << TmpV[j].Im(i) << ")    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, NDMAPSpinChainHamiltonian& H) 
{
  ComplexVector TmpV2 (H.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [H.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(H.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      ((AbstractHamiltonian*) &H)->LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chain->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chain->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Chain->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Chain->GetHilbertSpaceDimension() - 1][H.Chain->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

