////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of unconventional 3D toric code                  //
//                                                                            //
//                        last modification : 23/05/2017                      //
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


#ifndef UNCONVENTIONAL3DTORICCODEHAMILTONIAN_H
#define UNCONVENTIONAL3DTORICCODEHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Hamiltonian/AbstractHamiltonian.h"


#include <iostream>


using std::ostream;
class MathematicaOutput;


class Unconventional3DToricCodeHamiltonian : public AbstractHamiltonian
{

 protected:
  
  //pointer to Hilbert space of the associated system
  AbstractSpinChain* Chain;

  // total number of spins
  int NbrSpin;
  // number of spin chain along the x direction
  int NbrSpinX;
  // number of spin chain along the y direction
  int NbrSpinY;
  // number of spin chain along the z direction
  int NbrSpinZ;

  // true if periodic boundary conditions have to be used along the x direction
  bool PeriodicBoundaryConditionsX;
  // true if periodic boundary conditions have to be used along the x direction
  bool PeriodicBoundaryConditionsY;
  // true if periodic boundary conditions have to be used along the x direction
  bool PeriodicBoundaryConditionsZ;
  
  // global energy shift to apply to the system
  double HamiltonianShift;

 public:
   
   // default constructor
  //
  Unconventional3DToricCodeHamiltonian();

  // constructor from default data
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpinX = number of spin along the x direction
  // nbrSpinY = number of spin along the y direction
  // nbrSpinZ = number of spin along the z direction
  // periodicBoundaryConditionsX = true if periodic boundary conditions have to be usedalong the x direction
  // periodicBoundaryConditionsY = true if periodic boundary conditions have to be usedalong the y direction
  // periodicBoundaryConditionsZ = true if periodic boundary conditions have to be usedalong the z direction
  Unconventional3DToricCodeHamiltonian(AbstractSpinChain* chain, int nbrSpinX, int nbrSpinY, int nbrSpinZ,
				       bool periodicBoundaryConditionsX, bool periodicBoundaryConditionsY, bool periodicBoundaryConditionsZ);

  // destructor
  //
  ~Unconventional3DToricCodeHamiltonian();

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
 
  // apply the cube operator to one basis state
  //
  // cornerX = lower leftmost frontmost cube corner x coordinate
  // cornerY = lower leftmost frontmost cube corner y coordinate
  // cornerZ = lower leftmost frontmost cube corner z coordinate
  // vDestination = reference on the vector to which result has to be added
  // index = basis state index
  // coefficient = basis state coefficient
  virtual void ApplyCubeOperator(int cornerX, int cornerY, int cornerZ, RealVector& vDestination, int index, double coefficient);

  // apply the plaquette operator in the xy plane to one basis state
  //
  // cornerX = leftmost frontmost plaquette corner x coordinate
  // cornerY = leftmost frontmost plaquette corner y coordinate
  // cornerZ = plaquette z coordinate
  // vDestination = reference on the vector to which result has to be added
  // index = basis state index
  // coefficient = basis state coefficient
  virtual void ApplyPlaquetteOperator(int cornerX, int cornerY, int cornerZ, RealVector& vDestination, int index, double coefficient);

  // get a linearized position index from the 3d coordinates
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // zPosition = position along the z direction
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition, int zPosition);

};

// get a linearized position index from the 3d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// zPosition = position along the z direction
// return value = linearized index

inline int Unconventional3DToricCodeHamiltonian::GetLinearizedIndex(int xPosition, int yPosition, int zPosition)
{
  return (((xPosition * this->NbrSpinY) + yPosition) * this->NbrSpinZ) + zPosition;
}

#endif
