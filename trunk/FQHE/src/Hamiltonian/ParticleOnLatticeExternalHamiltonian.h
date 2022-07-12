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


#ifndef PARTICLEONLATTICEEXTERNALHAMILTONIAN_H
#define PARTICLEONLATTICEEXTERNALHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Tools/FQHESpectrum/LatticePhases.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnLatticeExternalHamiltonian : public AbstractQHEOnLatticeHamiltonian
{
 protected:

  // filename for external one-particle terms
  char *OneParticleTerms;
  
  // filename for external two-particle interactions
  char *TwoParticleTerms;

  // lattice structure (if known)
  LatticePhases *LatticeGeometry;
  
 public:

  // constructor for contact interactions on a square lattice
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrStates = number of quantum states
  // oneParticleTerms = file describing single particle terms
  // twoParticleTerms = file describing two-particle terms
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // hermitianFlag = flag indicating whether to use hermitian symmetry
  ParticleOnLatticeExternalHamiltonian(ParticleOnLattice* particles, int nbrParticles, int nbrStates, const char* oneParticleTerms, const char* twoParticleTerms, AbstractArchitecture* architecture, unsigned long memory = 0, char* precalculationFileName = 0, bool hoppingOnly = false, LatticePhases *latticeGeometry = NULL);

  // destructor
  //
  ~ParticleOnLatticeExternalHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  
  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnLatticeExternalHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeExternalHamiltonian& H);

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

};

#endif // PARTICLEONLATTICEEXTERNALHAMILTONIAN_H
