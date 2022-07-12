////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 19/10/2004                        //
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
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "Hamiltonian/PeriodicElectronHole3DHamiltonian.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"
#include "HilbertSpace/PeriodicThreeDOneParticle.h"

#include <iostream>
#include <math.h>

using std::ostream;
using std::cout;
using std::endl;

#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define COULOMBIAN_FACTOR 180.79

// each state is coded with 5 bits
#define NBRBITZ2 5
#define DELTAZ2  15

#define NBRBITY2 5
#define DELTAY2  15

#define NBRBITX2 5
#define DELTAX2  15

#define NBRBITZ1 5
#define DELTAZ1  15

#define NBRBITY1 5
#define DELTAY1  15

#define NBRBITX1 5
#define DELTAX1  15

#define HEX2 0x7fff 
#define HEX1 0x7fff
// end of 5 bit coding


#define NBRBIT2 (NBRBITX2 + NBRBITY2 + NBRBITZ2)
#define NBRBIT1 (NBRBITX1 + NBRBITY1 + NBRBITZ1)

#define SHIFTZ2 0
#define SHIFTY2 (SHIFTZ2 + NBRBITZ2)
#define SHIFTX2 (SHIFTY2 + NBRBITY2)
#define SHIFTZ1 (SHIFTX2 + NBRBITX2)
#define SHIFTY1 (SHIFTZ1 + NBRBITZ1)
#define SHIFTX1 (SHIFTY1 + NBRBITY1)
#define SHIFTZ1bis 0
#define SHIFTY1bis (SHIFTZ1bis + NBRBITZ1)
#define SHIFTX1bis (SHIFTY1bis + NBRBITZ2)

#define DELTAX ((DELTAX1 << SHIFTX1bis) | (DELTAY1 <<  SHIFTY1bis) | (DELTAZ1 << SHIFTZ1bis))
#define DELTAY ((DELTAX2 << SHIFTX2) | (DELTAY2 <<  SHIFTY2) | (DELTAZ2 << SHIFTZ2))

// constructor
//
// space = pointer to the Hilbert space of two particles
// Mex, Mey, Mez = effective masses in three directions of electron (in vacuum electron mass unit)
// Mhx, Mhy, Mhz = effective masses in three directions of hole (in vacuum electron mass unit)
// potentialElectron = pointer to the potential for electron
// potentialHole = pointer to the potential for hole
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
// dielectric = dielectric constant in the sample

PeriodicElectronHole3DHamiltonian::PeriodicElectronHole3DHamiltonian (PeriodicThreeDTwoParticles* space, double Mex, double Mey, double Mez, double Mhx, double Mhy, double Mhz, ThreeDConstantCellPotential* potentialElectron, ThreeDConstantCellPotential* potentialHole, double xSize, double ySize, double zSize, double dielectric)
{
  this->Space = space;
  PeriodicThreeDOneParticle* firstParticle = (PeriodicThreeDOneParticle*) this->Space->GetFirstParticleSpace ();
  PeriodicThreeDOneParticle* secondParticle = (PeriodicThreeDOneParticle*) this->Space->GetSecondParticleSpace ();
  this->NbrState1X = firstParticle->GetNbrStateX ();
  this->NbrState1Y = firstParticle->GetNbrStateY ();
  this->NbrState1Z = firstParticle->GetNbrStateZ ();
  this->NbrState2X = secondParticle->GetNbrStateX ();
  this->NbrState2Y = secondParticle->GetNbrStateY ();
  this->NbrState2Z = secondParticle->GetNbrStateZ ();

  if (((this->NbrState1X * 2 - 1) > (1 << NBRBITX1)) || ((this->NbrState1Y * 2 - 1) > (1 << NBRBITY1)) || ((this->NbrState1Z * 2 - 1) > (1 << NBRBITZ1)) || ((this->NbrState2X * 2 - 1) > (1 << NBRBITX2)) || ((this->NbrState2Y * 2 - 1) > (1 << NBRBITY2)) || ((this->NbrState2Z * 2 - 1) > (1 << NBRBITZ2)))
    {
      cout << "At least a number of states in a direction is too big.!" << endl;
      exit (1);
    }

  this->MakeConversionTable ();
  cout << "Evaluating the kinetic terms ..." << endl;
  this->EvaluateKineticTerm (Mex, Mey, Mez, Mhx, Mhy, Mhz, firstParticle, secondParticle, xSize, ySize, zSize);
  cout << "Evaluating the electron confinement terms ..." << endl;
  this->EvaluateConfinementTerm (potentialElectron, firstParticle, this->RealElectronConfinement, this->ImaginaryElectronConfinement);
  cout << "Evaluating the hole confinement terms ..." << endl;
  this->EvaluateConfinementTerm (potentialHole, secondParticle, this->RealHoleConfinement, this->ImaginaryHoleConfinement);
  cout << "Evaluating the Coulombian term ..." << endl;
  this->EvaluateCoulombianTerm (xSize, ySize, zSize, dielectric);
  firstParticle = NULL; secondParticle = NULL;
  cout << "End of the evaluation." << endl;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

PeriodicElectronHole3DHamiltonian::PeriodicElectronHole3DHamiltonian(const PeriodicElectronHole3DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->NbrState1X = NbrState1X;
  this->NbrState1Y = NbrState1Y;
  this->NbrState1Z = NbrState1Z;
  this->NbrState2X = NbrState2X;
  this->NbrState2Y = NbrState2Y;
  this->NbrState2Z = NbrState2Z;
  this->IToXY = hamiltonian.IToXY;
  this->XToI = hamiltonian.XToI;
  this->YToI = hamiltonian.YToI;
  this->KineticTerm = hamiltonian.KineticTerm;
  this->RealElectronConfinement = hamiltonian.RealElectronConfinement;
  this->ImaginaryElectronConfinement = hamiltonian.ImaginaryElectronConfinement;
  this->RealHoleConfinement = hamiltonian.RealHoleConfinement;
  this->ImaginaryHoleConfinement = hamiltonian.ImaginaryHoleConfinement;
  this->CoulombianTerm = hamiltonian.CoulombianTerm;
}

// destructor
//

PeriodicElectronHole3DHamiltonian::~ PeriodicElectronHole3DHamiltonian()
{
  //cout << "PeriodicElectronHole3DHamiltonian destructor is being called." << endl;
  delete   this->Space;
  //cout << "Destructor of Hilbert space is being called" << endl;
  delete[] this->IToXY;
  delete[] this->XToI;
  delete[] this->YToI;
  delete[] this->KineticTerm;
  delete[] this->RealElectronConfinement;
  delete[] this->ImaginaryElectronConfinement;
  delete[] this->RealHoleConfinement;
  delete[] this->ImaginaryHoleConfinement;
  delete[] this->CoulombianTerm;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* PeriodicElectronHole3DHamiltonian::Clone ()
{
  return new PeriodicElectronHole3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicElectronHole3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicElectronHole3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticTerm[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicElectronHole3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
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

Complex PeriodicElectronHole3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
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

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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

ComplexVector& PeriodicElectronHole3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{ 
  int dimension = this->Space->GetHilbertSpaceDimension ();
  int lastComponent = firstComponent + nbrComponent;
  int XY1 = 0, XY2 = 0; int Sum1 = 0;
  int X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;
  double TmpRe = 0.0, TmpIm = 0.0;
  int X1Shifted = 0, Y1Shifted = 0;

  for (int Index1 = firstComponent; Index1 < lastComponent; ++Index1)
    {
      XY1 = this->IToXY[Index1];
      Y1 = XY1 & HEX2; X1 = (XY1 >> NBRBIT2) & HEX1;
      Sum1 = X1 + Y1;
      X1Shifted = DELTAX - X1; Y1Shifted = DELTAY - Y1;      
      
      TmpRe = 0.0; TmpIm = 0.0;

      // kinetic term
      TmpRe += this->KineticTerm[Index1] * vSource.Re(Index1);
      TmpIm += this->KineticTerm[Index1] * vSource.Im(Index1);

      for (int Index2 = 0; Index2 < dimension; ++Index2)
	{
	  XY2 = this->IToXY[Index2];
	  Y2 = XY2 & HEX2; X2 = (XY2 >> NBRBIT2) & HEX1;
	  
	  double ReSource =  vSource.Re(Index2), ImSource = vSource.Im(Index2);
	  
	  // Coulombian term	  	  
	  if ((X2 + Y2) == Sum1)  
	    {	     
	      double tmpCoulomb = this->CoulombianTerm[this->XToI[X1Shifted + X2]];
	      TmpRe += tmpCoulomb * ReSource;
	      TmpIm += tmpCoulomb * ImSource;
	    }	  
	  
	  // confinement terms for the electrons
	  if (Y1 == Y2)
	    {
	      int tmpIndex = this->XToI[X1Shifted + X2];
	      double tmpERe = this->RealElectronConfinement[tmpIndex];
	      double tmpEIm = this->ImaginaryElectronConfinement[tmpIndex];
	      TmpRe += (tmpERe * ReSource - tmpEIm * ImSource);
	      TmpIm += (tmpERe * ImSource + tmpEIm * ReSource);
	    }
	  	  
	  // confinement terms for the holes
	  if (X1 == X2)
	    {
	      int tmpIndex = this->YToI[Y1Shifted + Y2];
	      double tmpERe = this->RealHoleConfinement[tmpIndex];
	      double tmpEIm = this->ImaginaryHoleConfinement[tmpIndex];	      
	      TmpRe += (tmpERe * ReSource - tmpEIm * ImSource);
	      TmpIm += (tmpERe * ImSource + tmpEIm * ReSource);
	    }	  
	}
      vDestination.Re(Index1) += TmpRe; vDestination.Im(Index1) += TmpIm;
    }

  return vDestination;
}
 
// determine the maximal value of the kenetic elements
//
// return = the wanted value

double PeriodicElectronHole3DHamiltonian::MaxKineticElement()
{
  double tmp = this->KineticTerm[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticTerm[i])
      tmp = this->KineticTerm[i];
  return tmp;
}

// make the conversion table to hexadecimal indices
//

void PeriodicElectronHole3DHamiltonian::MakeConversionTable ()
{
  int dimension = this->Space->GetHilbertSpaceDimension ();
  
  this->IToXY = new int [dimension];
  int index = 0;
  //cout << "I to X:" << endl;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    for (int n1 = 0; n1 < this->NbrState1Y; ++n1)      
      for (int p1 = 0; p1 < this->NbrState1Z; ++p1)	    	    
	for (int m2 = 0; m2 < this->NbrState2X; ++m2)
	  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)      
	    for (int p2 = 0; p2 < this->NbrState2Z; ++p2)	    
	      {
		this->IToXY[index] = (m1 << SHIFTX1) | (n1 << SHIFTY1) | (p1 << SHIFTZ1) | (m2 << SHIFTX2) | (n2 << SHIFTY2) | (p2 << SHIFTZ2);
		//cout << "Index: " << index << " " << m1 << " " << n1 << " " << p1 << " " << m2 << " " << n2 << " " << p2 << " " << this->IToXY[index] << endl;
		++index;
	      }
  
  int dimension1 = (1 << NBRBIT1);
  this->XToI = new int [dimension1];

  int LengthX1 = this->NbrState1X * 2 - 1; int LengthY1 = this->NbrState1Y * 2 - 1; int LengthZ1 = this->NbrState1Z * 2 - 1;
  int OriginX1 = this->NbrState1X - 1; int OriginY1 = this->NbrState1Y - 1; int OriginZ1 = this->NbrState1Z - 1; 

  int index1 = 0; int tmpX = 0; int tmpX1 = 0; int tmpY1 = 0;
  for (int i = 0; i < LengthX1; ++i)
    {
      tmpX1 = ((i - OriginX1 + DELTAX1) << SHIFTX1bis);
      for (int j = 0; j < LengthY1; ++j)
	{
	  tmpY1 = (tmpX1 | ((j - OriginY1 + DELTAY1) << SHIFTY1bis));
	  for (int k = 0; k < LengthZ1; ++k)
	    {
	      tmpX =  (tmpY1 | ((k - OriginZ1 + DELTAZ1) << SHIFTZ1bis));	     
	      this->XToI[tmpX] = index1;
	      ++index1;
	    }
	}
    }

  int dimension2 = (1 << NBRBIT2);
  this->YToI = new int [dimension2];

  int LengthX2 = this->NbrState2X * 2 - 1; int LengthY2 = this->NbrState2Y * 2 - 1; int LengthZ2 = this->NbrState2Z * 2 - 1;
  int OriginX2 = this->NbrState2X - 1; int OriginY2 = this->NbrState2Y - 1; int OriginZ2 = this->NbrState2Z - 1; 

  int index2 = 0; int tmpY = 0; int tmpX2 = 0; int tmpY2 = 0;
  for (int i = 0; i < LengthX2; ++i)
    {
      tmpX2 = ((i - OriginX2 + DELTAX2) << SHIFTX2);
      for (int j = 0; j < LengthY2; ++j)
	{
	  tmpY2 = (tmpX2 | ((j - OriginY2 + DELTAY2) << SHIFTY2));
	  for (int k = 0; k < LengthZ2; ++k)
	    {
	      tmpY =  (tmpY2 | ((k - OriginZ2 + DELTAZ2) << SHIFTZ2));	     
	      this->YToI[tmpY] = index2;
	      ++index2;
	    }
	}
    }

}

// evaluate the kinetic term
//
// mex, mey, mez = effective masses in three directions of electron (in vacuum electron mass unit)
// mhx, mhy, mhz = effective masses in three directions of hole (in vacuum electron mass unit)
// firstParticle = pointer to the first particle's Hilbert space
// secondParticle = pointer to the second particle's Hilbert space
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)

void PeriodicElectronHole3DHamiltonian::EvaluateKineticTerm (double mex, double mey, double mez, double mhx, double mhy, double mhz, PeriodicThreeDOneParticle* firstParticle, PeriodicThreeDOneParticle* secondParticle, double xSize, double ySize, double zSize)
{
  int dimension = this->Space->GetHilbertSpaceDimension ();
  this->KineticTerm = new double [dimension];
  int LowerImpulsion1X = firstParticle->GetLowerImpulsionX();
  int LowerImpulsion1Y = firstParticle->GetLowerImpulsionY();
  int LowerImpulsion1Z = firstParticle->GetLowerImpulsionZ();
  int LowerImpulsion2X = secondParticle->GetLowerImpulsionX();
  int LowerImpulsion2Y = secondParticle->GetLowerImpulsionY();
  int LowerImpulsion2Z = secondParticle->GetLowerImpulsionZ();
  
  double InvXFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mex * xSize * xSize);
  double InvYFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mey * ySize * ySize);
  double InvZFactor1 = PERIODIC_HAMILTONIAN_FACTOR / (mez * zSize * zSize);
  
  double InvXFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhx * xSize * xSize);
  double InvYFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhy * ySize * ySize);
  double InvZFactor2 = PERIODIC_HAMILTONIAN_FACTOR / (mhz * zSize * zSize);
  
  double* factorX1 = new double [this->NbrState1X]; double* factorY1 = new double [this->NbrState1Y]; double* factorZ1 = new double [this->NbrState1Z];
  double* factorX2 = new double [this->NbrState2X]; double* factorY2 = new double [this->NbrState2Y]; double* factorZ2 = new double [this->NbrState2Z];
  
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    factorX1[m1] = InvXFactor1 * ((double) (m1 + LowerImpulsion1X)) * ((double) (m1 + LowerImpulsion1X));
  for (int n1 = 0; n1 < this->NbrState1Y; ++n1)
    factorY1[n1] = InvYFactor1 * ((double) (n1 + LowerImpulsion1Y)) * ((double) (n1 + LowerImpulsion1Y));
  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
    factorZ1[p1] = InvZFactor1 * ((double) (p1 + LowerImpulsion1Z)) * ((double) (p1 + LowerImpulsion1Z));
  for (int m2 = 0; m2 < this->NbrState2X; ++m2)
    factorX2[m2] = InvXFactor2 * ((double) (m2 + LowerImpulsion2X)) * ((double) (m2 + LowerImpulsion2X));
  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
    factorY2[n2] = InvYFactor2 * ((double) (n2 + LowerImpulsion2Y)) * ((double) (n2 + LowerImpulsion2Y));
  for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
    factorZ2[p2] = InvZFactor2 * ((double) (p2 + LowerImpulsion2Z)) * ((double) (p2 + LowerImpulsion2Z));

  //cout << "Kinetic terms: " << endl;
  int index = 0; double FactorX1, FactorY1, FactorZ1, FactorX2, FactorY2;
  for (int m1 = 0; m1 < this->NbrState1X; ++m1)
    {
      FactorX1 = factorX1[m1];
      for (int n1 = 0; n1 < this->NbrState1Y; ++n1)	
	{
	  FactorY1 = FactorX1 + factorY1[n1];
	  for (int p1 = 0; p1 < this->NbrState1Z; ++p1)
	    {
	      FactorZ1 = FactorY1 + factorZ1[p1];	      
	      for (int m2 = 0; m2 < this->NbrState2X; ++m2)
		{
		  FactorX2 = FactorZ1 + factorX2[m2];
		  for (int n2 = 0; n2 < this->NbrState2Y; ++n2)
		    {
		      FactorY2 = FactorX2 + factorY2[n2];	  
		      for (int p2 = 0; p2 < this->NbrState2Z; ++p2)
			{
			  this->KineticTerm[index] = FactorY2 + factorZ2[p2];	 
			  //cout << index << " " <<  m1 << " " << n1 << " " << p1 << " " << m2 << " " << n2 << " " << p2 << " " << this->KineticTerm[index] << endl;
			  ++index;
			}
		    }
		}
	    }
	}
    }

  delete[] factorX1; delete[] factorY1; delete[] factorZ1;
  delete[] factorX2; delete[] factorY2; delete[] factorZ2;
}

// evaluate the confinement terms for electrons and holes
//
// potential = pointer to the potential for the considered carrier
// particle = pointer to the Hilbertspace for the considered carrier
// realConfinement = reference to 1D array of real elements of the wanted terms
// imaginaryConfinement = reference to 1D array of imaginary elements of the wanted terms

void PeriodicElectronHole3DHamiltonian::EvaluateConfinementTerm (ThreeDConstantCellPotential* potential, PeriodicThreeDOneParticle* particle, double* &realConfinement, double* &imaginaryConfinement)
{
  int NbrStateX = particle->GetNbrStateX (), NbrStateY = particle->GetNbrStateY (), NbrStateZ = particle->GetNbrStateZ ();
  int NbrCellX = potential->GetNbrCellX (), NbrCellY = potential->GetNbrCellY (), NbrCellZ = potential->GetNbrCellZ ();

  double** RealWaveFunctionOverlapX; double** RealWaveFunctionOverlapY; double** RealWaveFunctionOverlapZ; 
  double** ImaginaryWaveFunctionOverlapX; double** ImaginaryWaveFunctionOverlapY; double** ImaginaryWaveFunctionOverlapZ; 

  this->EvaluateWaveFunctionOverlap (NbrCellX, NbrStateX, RealWaveFunctionOverlapX, ImaginaryWaveFunctionOverlapX);
  this->EvaluateWaveFunctionOverlap (NbrCellY, NbrStateY, RealWaveFunctionOverlapY, ImaginaryWaveFunctionOverlapY);
  this->EvaluateWaveFunctionOverlap (NbrCellZ, NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ);

  int LengthX = (NbrStateX - 1) * 2 + 1; int LengthY = (NbrStateY - 1) * 2 + 1; int LengthZ = (NbrStateZ - 1) * 2 + 1;

  double*** TmpReal = new double** [LengthX];
  double*** TmpImaginary = new double** [LengthX];

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX;
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;

  for (int m = 0; m < LengthX; ++m)
    {
      TmpReal[m] = new double* [LengthY];
      TmpImaginary[m] = new double* [LengthY];
      TmpRealWaveFunctionOverlapX = RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [NbrCellZ];
	  TmpImaginary[m][n] = new double [NbrCellZ];
	  TmpRealWaveFunctionOverlapY = RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < NbrCellX; ++CellX)
		    {		      
		      TmpRe += potential->GetPotential (CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += potential->GetPotential (CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }
  
  int number = LengthX * LengthY * LengthZ;
  realConfinement = new double [number]; imaginaryConfinement = new double [number]; 
  
  int index = 0;
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    for (int n = 0; n < LengthY; ++n)
      {
	TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	for (int p = 0; p < LengthZ; ++p)
	  {
	    TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
	    TmpImaginaryWaveFunctionOverlapZ = ImaginaryWaveFunctionOverlapZ[p];
	    TmpRe = 0.0; TmpIm = 0.0;
	    for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
	      {
		TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
	      }

	    realConfinement[index] = TmpRe;
	    imaginaryConfinement[index] = TmpIm;   	    
	    ++index;
	  }
      }

  delete[] RealWaveFunctionOverlapX; delete[] RealWaveFunctionOverlapY; delete[] RealWaveFunctionOverlapZ;
  delete[] ImaginaryWaveFunctionOverlapX; delete[] ImaginaryWaveFunctionOverlapY; delete[] ImaginaryWaveFunctionOverlapZ;
  delete[] TmpReal; delete[] TmpImaginary;
}

// evaluate the Coulombian term
//
// xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
// dielectric = dielectric constant in the sample

void PeriodicElectronHole3DHamiltonian::EvaluateCoulombianTerm (double xSize, double ySize, double zSize, double dielectric)
{
  double squareX = xSize * xSize;
  double squareY = ySize * ySize;
  double squareZ = zSize * zSize;
  double volume = xSize * ySize * zSize;
  double factor = COULOMBIAN_FACTOR / (dielectric * volume * 4.0 * M_PI * M_PI);

  int lengthX = this->NbrState1X * 2 - 1, lengthY = this->NbrState1Y * 2 - 1, lengthZ = this->NbrState1Z * 2 - 1;
  int originX = this->NbrState1X - 1, originY = this->NbrState1Y - 1, originZ = this->NbrState1Z - 1;
  
  this->CoulombianTerm = new double [lengthX * lengthY * lengthZ];
  int index = 0;
  for (int i = 0; i < lengthX; ++i)
    for (int j = 0; j < lengthY; ++j)
      for (int k = 0; k < lengthZ; ++k)
	{
	  if ((i == originX) && (j == originY) && (k == originZ))
	    this->CoulombianTerm[index] = 0.0;
	  else
	    {
	      double tmp = ((double )(i - originX) * (i - originX)) / squareX + ((double) (j - originY) * (j - originY)) / squareY + ((double) (k - originZ) * (k - originZ)) / squareZ;
	      this->CoulombianTerm[index] = -factor / tmp;		
	    }
	  ++index;
	}
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool PeriodicElectronHole3DHamiltonian::EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Origin; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      Diff = 2.0 * M_PI * double (delta - Origin);
      Tmp = Diff / nbrStep;	
      Diff = 1.0 / Diff;	
      for (int i = 0; i < nbrStep; ++i)
	{
	  realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	  imaginaryArray[delta][i] = Diff * (-cos(Tmp * (i + 1)) + cos(Tmp * i));
	}    
    } 
  realArray[Origin] = new double [nbrStep];
  imaginaryArray[Origin] = new double [nbrStep];
  for (int i = 0; i < nbrStep; ++i)
    {
      realArray[Origin][i] = Tmp1;
      imaginaryArray[Origin][i] = 0.0;
    }
  for (int delta = Origin + 1; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];      
      for (int i = 0; i < nbrStep; ++i)
	{
	  realArray[delta][i] = realArray[Length - 1 - delta][i];
	  imaginaryArray[delta][i] = -imaginaryArray[Length - 1 - delta][i];
	}    	
    }
  return true;
}
