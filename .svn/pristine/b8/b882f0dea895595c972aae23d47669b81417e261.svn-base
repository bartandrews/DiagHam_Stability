////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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


#ifndef FERMIONONSPHEREDROPLET_H
#define FERMIONONSPHEREDROPLET_H


#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include <iostream>


class FermionOnSphere;
class RationalPolynomial;
class LongRationalPolynomial;


class FermionOnSphereDroplet :  public FermionOnSphere
{

  friend class FermionOnSphereHaldaneSymmetricBasis;
  friend class BosonOnSphereHaldaneBasisShort;

 protected:

  int NbrFluxes1;
  int MaxNbrHoles1;
  int MaxNbrParticles1;

  int NbrFluxes2;
  int MaxNbrHoles2;
  int MaxNbrParticles2;

 public:

  // default constructor
  //
  FermionOnSphereDroplet();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on twice the momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // nbrFluxes1 = condition for the number of fluxes in a droplet
  // maxNbrParticles1 = condition for the max number of particles in a droplet
  // maxNbrHoles1 = condition for the max number of holes in a droplet
  // nbrFluxes2 = secondary condition for the number of fluxes in a droplet
  // maxNbrParticles2 = secondary condition for the max number of particles in a droplet
  // maxNbrHoles2 = secondary condition for the max number of holes in a droplet
  // memory = amount of memory granted for precalculations
  FermionOnSphereDroplet (int nbrFermions, int totalLz, int lzMax, int nbrFluxes1, int maxNbrParticles1, int maxNbrHoles1, int nbrFluxes2, int maxNbrParticles2, int maxNbrHoles2, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereDroplet(const FermionOnSphereDroplet& fermions);

  // destructor
  //
  virtual ~FermionOnSphereDroplet ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereDroplet& operator = (const FermionOnSphereDroplet& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

};


#endif
