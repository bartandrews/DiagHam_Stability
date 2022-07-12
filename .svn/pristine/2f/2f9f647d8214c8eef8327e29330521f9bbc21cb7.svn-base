////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of two dimension spin model that could host             //
//                            a Read-Rezayi Z3 phase                          //
//                                                                            //
//                        last modification : 24/05/2018                      //
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


#ifndef TWODIMENSIONALRRHAMILTONIAN_H
#define TWODIMENSIONALRRHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/TwoDimensionalHeisenbergHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class TwoDimensionalRRHamiltonian : public TwoDimensionalHeisenbergHamiltonian
{

 protected:
  
  // amplitude of the Heisenberg coupling between nearest neighbors
  double J1Factor;
  // amplitude of the (S_i S_j)^2 nearest neighbor coupling
  double J2Factor;
  // amplitude of the (S_i S_j)^3 nearest neighbor coupling
  double J3Factor;
  // amplitude of the chiral term
  double JcFactor;
  // half the amplitude of the chiral term
  double HalfJcFactor;

  // assume periodic boundary conditions
  bool PeriodicBoundaryConditions;

 public:
   
   // default constructor
  //
  TwoDimensionalRRHamiltonian();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // j1Factor = amplitude of the Heisenberg coupling between nearest neighbors
  // j2Factor = amplitude of the (S_i S_j)^2 nearest neighbor coupling
  // j3Factor = amplitude of the (S_i S_j)^3 nearest neighbor coupling
  // jcFactor = amplitude of the chiral term
  // periodicBoundaryConditions = assume periodic boundary conditions
  TwoDimensionalRRHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, double j1Factor, 
			      double j2Factor, double j3Factor, double jcFactor, bool periodicBoundaryConditions);

  // destructor
  //
  ~TwoDimensionalRRHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

 protected:
 
  // evaluate all matrix elements
  //   
  virtual void EvaluateDiagonalMatrixElements();

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, RealVector& vDestination, double coefficient);

  // evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

  // evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
  //
  // i = linearized position of the first spin
  // j = linearized position of the second spin
  // k = linearized position of the second spin
  // index = index of the many-body state to act on
  // dimension = total Hilbert space dimension
  // vDestination = vector at which result has to be added
  // coefficient = global multiplicative coefficient
  virtual void EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient);

};

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

inline void TwoDimensionalRRHamiltonian::EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, RealVector& vDestination, double coefficient)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int pos = this->Chain->SmiSpj(i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      double TmpValue2 = 0.5 * TmpCoefficient * coefficient;
      vDestination[pos] += this->J1Factor * TmpValue2;
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      TmpValue2 *= TmpValue3;
      vDestination[pos] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
      TmpValue2 *= this->Chain->SziSzj(i, j, pos); 
      vDestination[pos] +=  this->J3Factor * TmpValue2; 
    }
  pos = this->Chain->SmiSpjSmiSpj(i, j, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      double TmpValue2 = 0.25 * TmpCoefficient * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;  
    }	  
  pos = this->Chain->SmiSpjSmiSpj(j, i, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      double TmpValue2 = 0.25 * TmpCoefficient * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;   
    }	  
  pos = this->Chain->SziSzjSmiSpj(i, j, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      double TmpValue2 = 0.5 * TmpCoefficient * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  0.5 * this->J3Factor * TmpValue2 * TmpCoefficient2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;
    }	  
}

// evaluate the off-diagonal contribution for one type of Hamiltonian terms ( all (S_i S_j)^n )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

inline void TwoDimensionalRRHamiltonian::EvaluateOffDiagonalPowerHeisenbergContribution(int i, int j, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  double TmpCoefficient2;
  int pos = this->Chain->SmiSpj(i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * coefficient;
      vDestination[pos] += this->J1Factor * TmpValue2;
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      TmpValue2 *= TmpValue3;
      vDestination[pos] += (this->J2Factor + (this->J3Factor * TmpValue3)) * TmpValue2;
      TmpValue2 *= this->Chain->SziSzj(i, j, pos); 
      vDestination[pos] +=  this->J3Factor * TmpValue2; 
    }
  pos = this->Chain->SmiSpjSmiSpj(i, j, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;  
    }	  
  pos = this->Chain->SmiSpjSmiSpj(j, i, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.25 * TmpCoefficient) * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      double TmpValue3 = this->Chain->SziSzj(i, j, index);
      vDestination[pos] +=  this->J3Factor * TmpValue2 * TmpValue3;
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;   
    }	  
  pos = this->Chain->SziSzjSmiSpj(i, j, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      Complex TmpValue2 = (0.5 * TmpCoefficient) * coefficient;
      vDestination[pos] += this->J2Factor * TmpValue2;	      
      int pos2 = this->Chain->SmiSpj(i, j, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      pos2 = this->Chain->SmiSpj(j, i, pos, TmpCoefficient2);
      if (pos2 != dimension)
	{
	  vDestination[pos2] +=  (0.5 * this->J3Factor * TmpCoefficient2) * TmpValue2;
	}
      TmpValue2 *= this->Chain->SziSzj(i, j, pos);
      vDestination[pos] +=  this->J3Factor * TmpValue2;
    }	  
}

// evaluate the off-diagonal chiral contribution for a single term ( S_i (S_j ^ S_k) )
//
// i = linearized position of the first spin
// j = linearized position of the second spin
// k = linearized position of the second spin
// index = index of the many-body state to act on
// dimension = total Hilbert space dimension
// vDestination = vector at which result has to be added
// coefficient = global multiplicative coefficient

inline void TwoDimensionalRRHamiltonian::EvaluateOffDiagonalChiralContribution(int i, int j, int k, int index, int dimension, ComplexVector& vDestination, Complex& coefficient)
{
  double TmpCoefficient;
  int pos = this->Chain->SpiSmjSzk(i, j, k, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(i, k, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, k, i, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, j, i, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(k, i, j, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re -= TmpCoefficient * coefficient.Im;
      vDestination[pos].Im += TmpCoefficient * coefficient.Re;
    }
  pos = this->Chain->SpiSmjSzk(j, i, k, index, TmpCoefficient);
  if (pos != dimension)
    {
      TmpCoefficient *= this->HalfJcFactor;
      vDestination[pos].Re += TmpCoefficient * coefficient.Im;
      vDestination[pos].Im -= TmpCoefficient * coefficient.Re;
    }
}

#endif
