////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 1 dimensions        //
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


#ifndef PERIODIC1DHAMILTONIAN_H
#define PERIODIC1DHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/PeriodicOneDOneParticle.h"
#include <iostream>


using std::ostream;


class MathematicaOutput;
class OneDConstantCellPotential;

class Periodic1DHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  PeriodicOneDOneParticle* Space;

  // wave function basis dimension 
  int NbrState;  

  int LowerImpulsion;   

  // system dimension (in Angstrom unit)
  double Size;

  // effective mass (in electron mass unit)
  double Mu;

  // cache for hamiltonian diagonal elements
  double* KineticElements;

  // hamiltonian elements containing only interaction terms
  double* RealHamiltonian;
  double* ImaginaryHamiltonian;

 public:

  // constructor from data
  //
  // space = Hilbert space
  // mu = effective mass
  // PotentialInput = pointer to a 1D potential with constant value in a cell
  // waveVector = wave vector of Bloch function
  Periodic1DHamiltonian(PeriodicOneDOneParticle* space, double mu, OneDConstantCellPotential* PotentialInput, double waveVector = 0.0);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  Periodic1DHamiltonian(const Periodic1DHamiltonian& hamiltonian);

  // destructor
  //
  ~Periodic1DHamiltonian();
  
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
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  
  Complex MatrixElement (RealVector& V1, RealVector& V2);


  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element

  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);
  
  // determine the maximal value of partial diagonal array
  //
  // return = the wanted value
  double MaxPartialDiagonalElement();

 private:

  // evaluate all interaction factors
  // 
  // potential = pointer to the potential
  // waveVector = wave vector of Bloch function
  void EvaluateInteractionFactors(OneDConstantCellPotential* &potential, double waveVector);

  // evaluate the plane wave function overlap
  //
  // potential = pointer to the potential
  // nbrState = number of states chosen for this direction
  // realArray = 2D array containing the real elements of the overlap
  // imaginaryArray = 2D array containing the imaginary elements of the overlap
  bool EvaluatePlaneWaveFunctionOverlap(OneDConstantCellPotential* &potential, int nbrState, double** &realArray, double** &imaginaryArray);
  
};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* Periodic1DHamiltonian::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int Periodic1DHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
