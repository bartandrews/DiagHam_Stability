////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of double triangle spin chain  hamiltonian           //
//                       with periodic boundary conditions                    //
//                                                                            //
//                        last modification : 19/06/2010                      //
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


#include "Hamiltonian/DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::endl;
using std::ostream;
using std::cout;


// constructor from default datas
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spins
// j1 = nearest neighbour coupling constant
// j2 = second nearest neighbour coupling constant
// djz1 = constant added to the nearest neighbour constant between spins along z
// djz2 = constant added to the second nearest neighbour constant between spins along z

DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j1, double j2, double djz1, double djz2)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->NbrTriangles = this->NbrSpin / 2;
  this->J1 = j1;
  this->HalfJ1 = this->J1 * 0.5;;
  this->J2 = j2;
  this->HalfJ2 = this->J2 * 0.5;;
  this->Jz1 = this->J1 + djz1;
  this->Jz2 = this->J2 + djz2;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::~DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (AbstractSpinChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian& DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::SetChain(AbstractSpinChainWithTranslations* chain)
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

int DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											      int firstComponent, int nbrComponent) 
{
  int dim = firstComponent + nbrComponent;
  int Dimension =  this->Chain->GetHilbertSpaceDimension();
  double Coef;
  int NbrTranslations;
  int pos;
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 2;
  Complex TmpZ;
  for (int i = firstComponent; i < dim; i++)
    {
      SpinPos = 0;
      Complex TmpCoefficient =  vSource[i];
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ1;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ1;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ1;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ1;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ2;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ2;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 3, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ2;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 3, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      Coef *= this->HalfJ2;
	      TmpZ = this->ExponentialTable[NbrTranslations];
	      TmpZ *= TmpCoefficient;
	      TmpZ *= Coef;
	      vDestination[pos] += TmpZ;
	    }
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ1;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ1;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, 0, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ1;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(0, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ1;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(SpinPos, 0, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ2;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(0, SpinPos, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ2;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ2;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      pos = this->Chain->SmiSpj(1, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  Coef *= this->HalfJ2;
	  TmpZ = this->ExponentialTable[NbrTranslations];
	  TmpZ *= TmpCoefficient;
	  TmpZ *= Coef;
	  vDestination[pos] += TmpZ;
	}
      vDestination[i] += this->SzSzContributions[i] * TmpCoefficient;
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors,
												      int firstComponent, int nbrComponent)
{
  int dim = firstComponent + nbrComponent;
  int Dimension =  this->Chain->GetHilbertSpaceDimension();
  double Coef;
  int NbrTranslations;
  int pos;
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 2;
  Complex* TmpCoefficients = new Complex [nbrVectors];
  for (int i = firstComponent; i < dim; i++)
    {
      SpinPos = 0;
      for (int j = 0; j < nbrVectors; ++j)
	TmpCoefficients[j] =  vSources[j][i];
      for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	{
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 2, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos, SpinPos + 2, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 2, SpinPos, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos + 3, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	    }
	  pos = this->Chain->SmiSpj(SpinPos + 3, SpinPos + 1, i, Coef, NbrTranslations);
	  if (pos != Dimension)
	    {
	      for (int j = 0; j <  nbrVectors; ++j)
		vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	    }
	}
      pos = this->Chain->SmiSpj(SpinPos, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, SpinPos, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, 0, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(0, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ1 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(SpinPos, 0, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(0, SpinPos, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(SpinPos + 1, 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	}
      pos = this->Chain->SmiSpj(1, SpinPos + 1, i, Coef, NbrTranslations);
      if (pos != Dimension)
	{
	  for (int j = 0; j <  nbrVectors; ++j)
	    vDestinations[j][pos] += this->HalfJ2 * Coef * TmpCoefficients[j];
	}
      for (int j = 0; j <  nbrVectors; ++j)
	vDestinations[j][i] += this->SzSzContributions[i] * TmpCoefficients[j];
    }
  delete[] TmpCoefficients;
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();

  // SzSz part
  int SpinPos = 0;
  int ReducedNbrSpin = this->NbrSpin - 2;
  for (int i = 0; i < dim; i++)
    {
       SpinPos = 0;
       // SzSz part
       double& Tmp = this->SzSzContributions[i];
       Tmp = 0.0;
       for (; SpinPos < ReducedNbrSpin; SpinPos += 2)
	 {
	   Tmp += this->Jz1 * (this->Chain->SziSzj(SpinPos, SpinPos + 1, i) + 
			       this->Chain->SziSzj(SpinPos + 1, SpinPos + 2, i));
	   Tmp += this->Jz2 * (this->Chain->SziSzj(SpinPos, SpinPos + 2, i) + 
			       this->Chain->SziSzj(SpinPos + 1, SpinPos + 3, i));
	 }
       this->SzSzContributions[i] += this->Jz1 * (this->Chain->SziSzj(SpinPos, SpinPos + 1, i)
						  + this->Chain->SziSzj(SpinPos + 1, 0, i));
       this->SzSzContributions[i] += this->Jz2 * (this->Chain->SziSzj(SpinPos, 0, i)
						  + this->Chain->SziSzj(SpinPos + 1, 1, i));      
    }
}

// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  this->ExponentialTable = new Complex [this->NbrSpin];
  double Coef = (2.0 * M_PI / (double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i].Re = this->CosinusTable[i];
      this->ExponentialTable[i].Im = this->SinusTable[i];
    }  
  // this->CosinusTable = new double [this->NbrTriangles];
//   this->SinusTable = new double [this->NbrTriangles];
//   this->ExponentialTable = new Complex [this->NbrTriangles];
//   double Coef = (2.0 * M_PI / (double) this->NbrTriangles) * ((double) this->Chain->GetMomentum());
//   for (int i = 0; i < this->NbrTriangles; ++i)
//     {
//       this->CosinusTable[i] = cos(Coef * ((double) i));
//       this->SinusTable[i] = sin(Coef * ((double) i));
//       this->ExponentialTable[i].Re = this->CosinusTable[i];
//       this->ExponentialTable[i].Im = this->SinusTable[i];
//     }
}
