////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk using the Haldane basis            //
//                                                                            //
//                        last modification : 01/07/2008                      //
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


#ifndef FERMIONONDISKHALDANEBASIS_H
#define FERMIONONDISKHALDANEBASIS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include <iostream>


class FermionOnDisk;


class FermionOnDiskHaldaneBasis :  public FermionOnSphereHaldaneBasis
{
  
 protected:

 public:

  // default constructor
  //
  FermionOnDiskHaldaneBasis();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on the angular momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax =  the maximum Lz value reached by a fermion
  // referenceState = array that describes the reference state to start from
  // memory = amount of memory granted for precalculations
  FermionOnDiskHaldaneBasis (int nbrFermions, int& totalLz, int lzMax, int* referenceState, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnDiskHaldaneBasis (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnDiskHaldaneBasis(const FermionOnDiskHaldaneBasis& fermions);

  // destructor
  //
  virtual ~FermionOnDiskHaldaneBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnDiskHaldaneBasis& operator = (const FermionOnDiskHaldaneBasis& fermions);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
  //
  // state = reference to the state to convert
  // shift = shift to apply to each orbitals
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // radius = radius of the A disk
  // groundState = reference on the total system ground state
  // shift = shift to apply to each orbitals
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)  
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixRealSpacePartition (int nbrFermionSector, int lzSector,  double radius , RealVector& groundState, int shift = 0, AbstractArchitecture* architecture = 0);
   
	
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
  // and computed from precalculated particle entanglement matrix
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // radius = radius of the A disk
  // shift = shift to apply to each orbitals
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  RealMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, double radius, int shift, RealMatrix& entanglementMatrix);
  
 protected:

  // core part of the convertion of a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // sqrtCoefficients = array that contains the normalization coefficients
  // invSqrtCoefficients = array that contains the inverts of the normalization coefficients
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomialCore(RealVector& state, double* sqrtCoefficients, double* invSqrtCoefficients, long reference, bool symmetryFactor);

};

#endif


