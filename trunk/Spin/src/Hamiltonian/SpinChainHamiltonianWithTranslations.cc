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


#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "Architecture/ArchitectureOperation/GenericHamiltonianPrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

SpinChainHamiltonianWithTranslations::SpinChainHamiltonianWithTranslations()
{
  this->PrecalculationShift = 0;
  this->FastMultiplicationFlag = false;
  this->FastMultiplicationStep = 1;
  this->NbrInteractionPerComponent = 0;
  this->NbrBalancedTasks = 0;
  this->LoadBalancingArray = 0;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->HermitianSymmetryFlag = false;
  this->Architecture = 0;
}

// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// j = coupling constant between spin
// nnCoupling = term to add to ZZ nearest-neighbour interaction
// nnnCoupling = nearest-neighbour interaction in Z direction

SpinChainHamiltonianWithTranslations::SpinChainHamiltonianWithTranslations(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double nnCoupling, double nnnCoupling)
{
  //cout << "SpinChainHamiltonianWithTranslations" << endl;
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->J = j;
  this->HalfJ = this->J * 0.5;
  this->Jz = this->J;
  if (nnCoupling != 0)
    {
      this->Jz += nnCoupling;
      cout << "Adjusting ZZ interaction: " << this->Jz << endl;
    }

  this->NNNCoupling = nnnCoupling;
  if (this->NNNCoupling != 0)
    cout << "Adding NNN interaction: " << this->NNNCoupling << endl;
  
  this->HermitianSymmetryFlag = false;
  this->Architecture = 0;
  this->FastMultiplicationFlag = false;
  this->FastMultiplicationStep = 1;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainHamiltonianWithTranslations::~SpinChainHamiltonianWithTranslations() 
{
  delete[] this->SzSzContributions;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  delete[] this->ExponentialTable;
  if (this->FastMultiplicationFlag == true)
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
	  if (this->NbrInteractionPerComponent[i] > 0)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	}
      delete[] this->InteractionPerComponentCoefficient;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainHamiltonianWithTranslations::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainHamiltonianWithTranslations::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

SpinChainHamiltonianWithTranslations& SpinChainHamiltonianWithTranslations::SetChain(AbstractSpinChainWithTranslations* chain)
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

int SpinChainHamiltonianWithTranslations::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainHamiltonianWithTranslations::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool SpinChainHamiltonianWithTranslations::SavePrecalculation (char* fileName)
{
  return false;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SpinChainHamiltonianWithTranslations::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex(0);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SpinChainHamiltonianWithTranslations::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
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
	  TmpZ.Re = Coef * ((V2.Re(i) * this->CosinusTable[NbrTranslation]) -
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
  return Z;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& SpinChainHamiltonianWithTranslations::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
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
  return vDestination;
}
 
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* SpinChainHamiltonianWithTranslations::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									      int firstComponent, int nbrComponent)
{
  double Coef;
  int NbrTranslation;
  int pos;
  int MaxPos = this->NbrSpin - 1;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; i++)
    {
      for (int k = 0; k < nbrVectors; ++k)
	{
	  vDestinations[k].Re(i) += this->SzSzContributions[i] * vSources[k].Re(i);
	  vDestinations[k].Im(i) += this->SzSzContributions[i] * vSources[k].Im(i);
	}
      for (int j = 0; j < MaxPos; j++)
	{
	  pos = this->Chain->SmiSpj(j, j + 1, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k].Re(pos) += Coef * ((vSources[k].Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSources[k].Im(i) * this->SinusTable[NbrTranslation]));
		  vDestinations[k].Im(pos) += Coef * ((vSources[k].Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSources[k].Im(i) * this->CosinusTable[NbrTranslation]));
		}
	    }
	  pos = this->Chain->SmiSpj(j + 1, j, i, Coef, NbrTranslation);
	  if (pos != this->Chain->GetHilbertSpaceDimension())
	    {
	      Coef *= this->HalfJ;
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k].Re(pos) += Coef * ((vSources[k].Re(i) * this->CosinusTable[NbrTranslation]) -
						      (vSources[k].Im(i) * this->SinusTable[NbrTranslation]));
		  vDestinations[k].Im(pos) += Coef * ((vSources[k].Re(i) * this->SinusTable[NbrTranslation]) +
						      (vSources[k].Im(i) * this->CosinusTable[NbrTranslation]));
		}
	    }
	}    
      pos = this->Chain->SmiSpj(MaxPos, 0, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k].Re(pos) += Coef * ((vSources[k].Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSources[k].Im(i) * this->SinusTable[NbrTranslation]));
	      vDestinations[k].Im(pos) += Coef * ((vSources[k].Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSources[k].Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	}
      pos = this->Chain->SmiSpj(0, MaxPos, i, Coef, NbrTranslation);
      if (pos != this->Chain->GetHilbertSpaceDimension())
	{
	  Coef *= this->HalfJ;
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k].Re(pos) += Coef * ((vSources[k].Re(i) * this->CosinusTable[NbrTranslation]) -
						  (vSources[k].Im(i) * this->SinusTable[NbrTranslation]));
	      vDestinations[k].Im(pos) += Coef * ((vSources[k].Re(i) * this->SinusTable[NbrTranslation]) +
						  (vSources[k].Im(i) * this->CosinusTable[NbrTranslation]));
	    }
	}      
    }
  return vDestinations;
}


// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> SpinChainHamiltonianWithTranslations::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void SpinChainHamiltonianWithTranslations::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  this->ExponentialTable = new Complex [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i].Re = this->CosinusTable[i];
      this->ExponentialTable[i].Im = this->SinusTable[i];
    }
}

// evaluate diagonal matrix elements
// 

void SpinChainHamiltonianWithTranslations::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  for (int i = 0; i < dim; i++)
    {
      this->SzSzContributions[i] = 0.0;
      for (int j = 0; j < (this->NbrSpin - 1); j++)
	{
	  this->SzSzContributions[i] += this->Chain->SziSzj(j, j + 1, i);
	}
      this->SzSzContributions[i] += this->Chain->SziSzj(this->NbrSpin - 1, 0, i);
      this->SzSzContributions[i] *= this->Jz;
      if (this->NNNCoupling != 0)
        {
          for (int j = 0; j < (this->NbrSpin - 2); j++)
   	     this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(j, j + 2, i);
   	  this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(this->NbrSpin - 2, 0, i);
   	  this->SzSzContributions[i] += this->NNNCoupling * this->Chain->SziSzj(this->NbrSpin - 1, 1, i);
        }
    }
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long SpinChainHamiltonianWithTranslations::FastMultiplicationMemory(long allowedMemory)
{
  if (this->Architecture == 0)
    {
      cout << "error, your architecture was not set. You cannot use fast-multiplication" << endl;
      return 0l;
    }
  else
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
      
      GenericHamiltonianPrecalculationOperation Operation(this);
      Operation.ApplyOperation(this->Architecture);
 
      if (this->Architecture->GetOptimizedTypicalRange(this->NbrInteractionPerComponent, MinIndex, MaxIndex) == true)
	{
	  this->PrecalculationShift = (int) MinIndex;
	  EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
	  cout << "distributed calculations have been reoptimized" << endl;
	}  
      if (allowedMemory == 0l)
	{
	  delete[] this->NbrInteractionPerComponent;
	  this->NbrInteractionPerComponent = 0;
	  return 0l;
	}
      long Memory = 0;
      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
	{
	  Memory += this->NbrInteractionPerComponent[i];
	}
      cout << "nbr interaction = " << Memory << endl;
      long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension;
      if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	{
	  this->FastMultiplicationStep = 1;
	  int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	    {
	      ++this->FastMultiplicationStep;
	      ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	      if (this->Chain->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
		++ReducedSpaceDimension;
	      TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension;
	      Memory = 0;
	      for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
		Memory += this->NbrInteractionPerComponent[i];
	    }
	  Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
	  long ResidualMemory = allowedMemory - Memory;
	  if (ResidualMemory > 0)
	    {
	      int TotalReducedSpaceDimension = ReducedSpaceDimension;
	      int* TmpNbrInteractionPerComponent = new int [TotalReducedSpaceDimension];
	      int i = 0;
	      int Pos = 0;
	      for (; i < ReducedSpaceDimension; ++i)
		{
		  TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[Pos];
		  Pos += this->FastMultiplicationStep;
		}
	      delete[] this->NbrInteractionPerComponent;
	      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	    }
	  else
	    {
	      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
	      for (int i = 0; i < ReducedSpaceDimension; ++i)
		TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
	      delete[] this->NbrInteractionPerComponent;
	      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	    }
	}
      else
	{
	  Memory = ((sizeof (int*) + sizeof (int) + sizeof(Complex*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
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
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long SpinChainHamiltonianWithTranslations::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0l;
  AbstractSpinChainWithTranslations* TmpSpinChain = (AbstractSpinChainWithTranslations*) this->Chain->Clone();
  int LastComponent = lastComponent + firstComponent;
  this->EvaluateFastMultiplicationMemoryComponent(TmpSpinChain, firstComponent, LastComponent, Memory);
  delete TmpSpinChain;
  return Memory;
}

// enable fast multiplication algorithm
//

void SpinChainHamiltonianWithTranslations::EnableFastMultiplication()
{
  if (this->Architecture == 0)
    {
      cout << "error, your architecture was not set. You cannot use fast-multiplication" << endl;
    }
  else
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int* TmpIndexArray;
      Complex* TmpCoefficientArray;
      long Pos;
      timeval TotalStartingTime2;
      timeval TotalEndingTime2;
      double Dt2;
      gettimeofday (&(TotalStartingTime2), 0);
      int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedSpaceDimension;
      this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
      this->InteractionPerComponentCoefficient = new Complex* [ReducedSpaceDimension];
      
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	{
	  this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
	  this->InteractionPerComponentCoefficient[i] = new Complex [this->NbrInteractionPerComponent[i]];
	}
      
      GenericHamiltonianPrecalculationOperation Operation(this, false);
      Operation.ApplyOperation(this->Architecture);
      
      this->FastMultiplicationFlag = true;
      gettimeofday (&(TotalEndingTime2), 0);
      cout << "------------------------------------------------------------------" << endl << endl;;
      Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
	((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
      cout << "time = " << Dt2 << endl;
    }
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void SpinChainHamiltonianWithTranslations::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  AbstractSpinChainWithTranslations* TmpSpinChain = (AbstractSpinChainWithTranslations*) this->Chain->Clone();

  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      long TotalPos = 0;
      this->EvaluateFastMultiplicationComponent(TmpSpinChain, i, this->InteractionPerComponentIndex[Pos], 
						this->InteractionPerComponentCoefficient[Pos], TotalPos);
      ++Pos;
    }
  delete TmpSpinChain;
}

// core part of the FastMultiplication method
// 
// chain = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray  

void SpinChainHamiltonianWithTranslations::EvaluateFastMultiplicationComponent(AbstractSpinChain* chain, int index, 
									       int* indexArray, Complex* coefficientArray, long& position)
{
  cout << "using dummy SpinChainHamiltonianWithTranslations::EvaluateFastMultiplicationComponent" << endl;
}

// core part of the PartialFastMultiplicationMemory
// 
// chain = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations  

void SpinChainHamiltonianWithTranslations::EvaluateFastMultiplicationMemoryComponent(AbstractSpinChain* chain, int firstComponent, int lastComponent, long& memory)
{
  cout << "using dummy SpinChainHamiltonianWithTranslations::EvaluateFastMultiplicationMemoryComponent" << endl;
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
//
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
// return value = true if no error occured

bool SpinChainHamiltonianWithTranslations::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  long MinIndex;
  long MaxIndex;
  if (this->Architecture == 0)
    {
       return false;
    } 
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  if ((this->NbrInteractionPerComponent != 0) && (this->FastMultiplicationStep != 0))
    {
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;

      if ((this->LoadBalancingArray == 0) || (this->NbrBalancedTasks != nbrTasks))
	{
	  if (this->LoadBalancingArray != 0)
	    delete [] this->LoadBalancingArray;
	  long* SegmentSize = new long[nbrTasks];
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  this->NbrBalancedTasks = nbrTasks;
	  long TmpNbrElement = 0l;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    TmpNbrElement += this->NbrInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement / ((long) nbrTasks);
	  TmpNbrElement = 0l;
	  int Pos = 0;
	  this->LoadBalancingArray[0] = MinIndex;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement += NbrInteractionPerComponent[i];
	      if (TmpNbrElement > TmpNbrPerSegment)
		{
		  SegmentSize[Pos] = TmpNbrElement;
		  this->LoadBalancingArray[Pos + 1] = MinIndex + (i * this->FastMultiplicationStep);
		  TmpNbrElement = 0l;
		  ++Pos;
		}
	    }
	  while (Pos < (nbrTasks - 1))
	    {
	      LoadBalancingArray[Pos + 1] = MaxIndex + 1;
	      SegmentSize[Pos] = 0;
	      ++Pos;
	    }
	  this->LoadBalancingArray[nbrTasks] = MaxIndex + 1;
	  SegmentSize[nbrTasks - 1] = TmpNbrElement;
	  
	  cout << "LoadBalancingArray=[ (" << (this->LoadBalancingArray[1] - this->LoadBalancingArray[0]) << ", " << SegmentSize[0] <<")";
	  for (int i = 1; i < nbrTasks; ++i)
	    cout << " (" << (this->LoadBalancingArray[i + 1] - this->LoadBalancingArray[i]) << ", " << SegmentSize[i] << ")";
	  cout << "]"<< endl;
	  delete[] SegmentSize;
	}
    }
  else
    {
      return false;
    }
  segmentIndices = this->LoadBalancingArray;
  return true;
}

