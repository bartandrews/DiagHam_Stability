////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//   to particles with contact interactions on a lattice in magnetic field    //
//                                                                            //
//                      last modification : 13/02/2008                        //
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


#ifndef PARTICLEONLATTICEDELTAHAMILTONIAN_H
#define PARTICLEONLATTICEDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnLatticeDeltaHamiltonian : public AbstractQHEOnLatticeHamiltonian
{
 protected:
  // strength of on-site delta-interaction
  double ContactInteractionU;

  // strength of delta-potential at the origin
  double DeltaPotential;

  // strength of random potential at individual sites
  double RandomPotential;

  // flag for reversed hopping
  bool ReverseHopping;

  // flag for cylinder boundary condition
  bool CylinderGeometry;


 public:

  // constructor for contact interactions on a square lattice
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lx = length of simulation cell in x-direction
  // ly = length of simulation cell in y-direction
  // nbrFluxQuanta = number of flux quanta piercing the simulation cell
  // contactInteractionU = strength of on-site delta interaction
  // reverseHopping = flag to indicate if sign of hopping terms should be reversed
  // deltaPotential = strength of a delta potential at site (0,0)
  // randomPotential = strength of a random on-site potential
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // hermitianFlag = flag indicating whether to use hermitian symmetry
  // cylinder_geometry = flag indicating whether to omit periodic boundary condition in the x-direction
  ParticleOnLatticeDeltaHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, double deltaPotential, double randomPotential, AbstractArchitecture* architecture, int nbrBody = 2, unsigned long memory = 0, char* precalculationFileName = 0, bool hermitianFlag = false, bool cylinder_geometry=false);

  // destructor
  //
  ~ParticleOnLatticeDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  
  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnLatticeDeltaHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeDeltaHamiltonian& H);

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

};

#endif // PARTICLEONLATTICEDELTAHAMILTONIAN_H
