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


#ifndef PERIODICELECTRONHOLE3DHAMILTONIAN_H
#define PERIODICELECTRONHOLE3DHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/PeriodicThreeDTwoParticles.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class ThreeDConstantCellPotential;

class PeriodicElectronHole3DHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  PeriodicThreeDTwoParticles* Space;

  // number of states in each direction for the first particle
  int NbrState1X;
  int NbrState1Y;
  int NbrState1Z;
  // number of states in each direction for the second particle
  int NbrState2X;
  int NbrState2Y;
  int NbrState2Z;

  // table of conversions to hexadecimal indices
  int* IToXY;
  int* XToI;
  int* YToI;
  
  double* KineticTerm;
  double* RealElectronConfinement;
  double* ImaginaryElectronConfinement;
  double* RealHoleConfinement;
  double* ImaginaryHoleConfinement;
  double* CoulombianTerm;

 public:

  // constructor
  //
  // space = pointer to the Hilbert space of two particles
  // Mex, Mey, Mez = effective masses in three directions of electron (in vacuum electron mass unit)
  // Mhx, Mhy, Mhz = effective masses in three directions of hole (in vacuum electron mass unit)
  // potentialElectron = pointer to the potential for electron
  // potentialHole = pointer to the potential for hole
  // xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
  // dielectric = dielectric constant in the sample
  PeriodicElectronHole3DHamiltonian (PeriodicThreeDTwoParticles* space, double Mex, double Mey, double Mez, double Mhx, double Mhy, double Mhz, ThreeDConstantCellPotential* potentialElectron, ThreeDConstantCellPotential* potentialHole, double xSize, double ySize, double zSize, double dielectric);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  PeriodicElectronHole3DHamiltonian(const PeriodicElectronHole3DHamiltonian& hamiltonian);

  // destructor
  //
  ~PeriodicElectronHole3DHamiltonian();
  
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
  
  // determine the maximal value of the kenetic elements
  //
  // return = the wanted value
  double MaxKineticElement();

 private:

  // make the conversion table to hexadecimal indices
  //
  void MakeConversionTable ();

  // evaluate the kinetic term
  //
  // mex, mey, mez = effective masses in three directions of electron (in vacuum electron mass unit)
  // mhx, mhy, mhz = effective masses in three directions of hole (in vacuum electron mass unit)
  // firstParticle = pointer to the first particle's Hilbert space
  // secondParticle = pointer to the second particle's Hilbert space
  // xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
  void EvaluateKineticTerm (double mex, double mey, double mez, double mhx, double mhy, double mhz, PeriodicThreeDOneParticle* firstParticle, PeriodicThreeDOneParticle* secondParticle, double xSize, double ySize, double zSize);

  // evaluate the confinement terms for electrons and holes
  //
  // potential = pointer to the potential for the considered carrier
  // particle = pointer to the Hilbertspace for the considered carrier
  // realConfinement = reference to 1D array of real elements of the wanted terms
  // imaginaryConfinement = reference to 1D array of imaginary elements of the wanted terms
  void EvaluateConfinementTerm (ThreeDConstantCellPotential* potential, PeriodicThreeDOneParticle* particle, double* &realConfinement, double* &imaginaryConfinement);

  // evaluate the Coulombian term
  //
  // xSize, ySize, zSize = sizes of the sample in three direction (in Angstrom unit)
  // dielectric = dielectric constant in the sample
  void EvaluateCoulombianTerm (double xSize, double ySize, double zSize, double dielectric);

  // evaluate the wave function overlap
  //
  // nbrStep = number of steps in the given direction
  // nbrState = number of states chosen for this direction
  // realArray = 2D array containing the real elements of the overlap
  // imaginaryArray = 2D array containing the imaginary elements of the overlap
  bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* PeriodicElectronHole3DHamiltonian::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int PeriodicElectronHole3DHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
