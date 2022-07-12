////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 24/11/2003                        //
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
#include "Hamiltonian/Periodic1DHamiltonian.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/Potential/OneDConstantCellPotential.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;


#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define BLOCH_FACTOR 7.644


// constructor from data
//
// space = Hilbert space
// mu = effective mass
// PotentialInput = pointer to a 1D potential with constant value in a cell
// waveVector = wave vector of Bloch function

Periodic1DHamiltonian::Periodic1DHamiltonian(PeriodicOneDOneParticle* space, double mu, OneDConstantCellPotential* PotentialInput, double waveVector)
{
  this->Space = space;
  this->NbrState = this->Space->GetNbrState();
  this->LowerImpulsion = this->Space->GetLowerImpulsion();
  this->Size = PotentialInput->GetSize();
  this->Mu = mu;  
  cout << "Evaluating confinement potential ..." << endl;
  this->EvaluateInteractionFactors(PotentialInput, waveVector);
  cout << "End of confinement potential evaluation." << endl;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

Periodic1DHamiltonian::Periodic1DHamiltonian(const Periodic1DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->Size = hamiltonian.Size;
  this->Mu = hamiltonian.Mu;
  this->NbrState = this->Space->GetNbrState();
  this->LowerImpulsion = this->Space->GetLowerImpulsion();
  this->KineticElements = hamiltonian.KineticElements;
  this->RealHamiltonian = hamiltonian.RealHamiltonian;
  this->ImaginaryHamiltonian = this->ImaginaryHamiltonian;
}

// destructor
//

Periodic1DHamiltonian::~ Periodic1DHamiltonian()
{
  delete[] this->KineticElements;
  delete[] this->RealHamiltonian;
  delete[] this->ImaginaryHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* Periodic1DHamiltonian::Clone ()
{
  return new Periodic1DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void Periodic1DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Periodic1DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex Periodic1DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  int dim = this->Space->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex Periodic1DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& Periodic1DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& Periodic1DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

ComplexVector& Periodic1DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& Periodic1DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  int lastComponent = firstComponent + nbrComponent;
  int Origin = this->NbrState - 1; 
  int m1, m2;
  int Index;

  double TmpRe = 0.0; double TmpIm = 0.0;
  int Limit = 0; int Length = (this->NbrState - 1) * 2 + 1;   
  
  for (m1 = firstComponent; m1 < lastComponent; ++m1)
    {
      vDestination.Re(m1) += vSource.Re(m1) * this->KineticElements[m1];
      vDestination.Im(m1) += vSource.Im(m1) * this->KineticElements[m1];   
      Limit = Length - m1;
      m2 = 0;
      TmpRe = 0.0; TmpIm = 0.0;
      for (Index = -m1 + Origin; Index < Limit; ++Index)
	{
	  TmpRe += (vSource.Re(m2) * this->RealHamiltonian[Index] - vSource.Im(m2) * this->ImaginaryHamiltonian[Index]);
	  TmpIm += (vSource.Re(m2) * this->ImaginaryHamiltonian[Index] + vSource.Im(m2) * this->RealHamiltonian[Index]);  	  
	  ++m2;
	}
      vDestination.Re(m1) += TmpRe; vDestination.Im(m1) += TmpIm;
    }
  return vDestination;
}

// evaluate all interaction factors
// 
// potential = pointer to the potential
// waveVector = wave vector of Bloch function
void Periodic1DHamiltonian::EvaluateInteractionFactors(OneDConstantCellPotential* &potential, double waveVector)
{
  cout << "size = " << this->Size << endl;
  double** RealWaveFunctionOverlap; double** ImaginaryWaveFunctionOverlap;
  if (!this->EvaluatePlaneWaveFunctionOverlap(potential, this->NbrState, RealWaveFunctionOverlap, ImaginaryWaveFunctionOverlap))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;
  double InvFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mu * this->Size * this->Size);

  double ShiftSquareK = BLOCH_FACTOR * (waveVector * waveVector / (2.0 * this->Mu));
  double ShiftK = BLOCH_FACTOR * waveVector * 2.0 * M_PI/ (this->Mu * this->Size); 
  this->KineticElements = new double[this->Space->GetHilbertSpaceDimension ()];

  for (int i = 0; i < this->NbrState; ++i)
    {	      
      this->KineticElements[i] = double((i + this->LowerImpulsion) * (i + this->LowerImpulsion)) * InvFactor;	  
      this->KineticElements[i] += (ShiftSquareK + ShiftK * double(i + this->LowerImpulsion));
      ++i;
    }

  int Length = (this->NbrState - 1) * 2 + 1; 

  this->RealHamiltonian = new double [Length];
  this->ImaginaryHamiltonian = new double [Length];
  double* TmpRealWaveFunctionOverlap;
  double* TmpImaginaryWaveFunctionOverlap;
  int nbrCell = potential->GetNumberCell();
  double TmpRe = 0.0, TmpIm = 0.0;
  for (int m = 0; m < Length; ++m)
    {
      TmpRealWaveFunctionOverlap = RealWaveFunctionOverlap[m];
      TmpImaginaryWaveFunctionOverlap = ImaginaryWaveFunctionOverlap[m];
      TmpRe = 0.0; TmpIm = 0.0;
      for (int i = 0; i < nbrCell; ++i)
	{
	  TmpRe += (TmpRealWaveFunctionOverlap[i] * potential->GetPotential(i));
	  TmpIm += (TmpImaginaryWaveFunctionOverlap[i] * potential->GetPotential(i));
	}
      this->RealHamiltonian[m] = TmpRe; this->ImaginaryHamiltonian[m] = TmpIm;
    }
  delete[] RealWaveFunctionOverlap; delete[] ImaginaryWaveFunctionOverlap;
}

// evaluate the plane wave function overlap
//
// potential = pointer to the potential
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool Periodic1DHamiltonian::EvaluatePlaneWaveFunctionOverlap(OneDConstantCellPotential* &potential, int nbrState, double** &realArray, double** &imaginaryArray)
{
  int Length = (nbrState - 1) * 2 + 1;
  int Origin = nbrState - 1;
  int nbr = potential->GetNumberCell();
  double* Position = new double [nbr + 1];
  Position[0] = 0.0;
  for (int k = 0; k < nbr; ++k)    
    Position[k + 1] = Position[k] + potential->GetCellWidth(k);    

  realArray = new double* [Length];
  imaginaryArray = new double* [Length];

  double Diff = 0.0, Tmp = 0.0; int delta = 0;
  for (int m = 0; m < Length; ++m)
    {
      delta = m - Origin;
      realArray[m] = new double [nbr];
      imaginaryArray[m] = new double [nbr];
      if (delta != 0)
	{
	  Diff = 2.0 * M_PI * double(delta);
	  Tmp = Diff / this->Size;
	  Diff = 1.0 / Diff;
	  for (int k = 0; k < nbr; ++k)
	    {
	      realArray[m][k] = Diff * (sin(Tmp * Position[k + 1]) - sin(Tmp * Position[k]));
	      imaginaryArray[m][k] = Diff * (cos(Tmp * Position[k]) - cos(Tmp * Position[k + 1]));
	    }
	}
      else 
	{
	  for (int k = 0; k < nbr; ++k)
	    {
	      realArray[Origin][k] = (Position[k + 1] - Position[k]) / this->Size;
	      imaginaryArray[Origin][k] = 0.0;
	    }
	}
    }
  return true;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double Periodic1DHamiltonian::MaxPartialDiagonalElement()
{
  double tmp = this->KineticElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticElements[i])
      tmp = this->KineticElements[i];
  return tmp;
}
