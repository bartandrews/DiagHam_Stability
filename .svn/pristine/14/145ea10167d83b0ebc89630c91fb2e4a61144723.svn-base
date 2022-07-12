////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 23/07/2010                      //
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


#ifndef PARTICLEONTORUSTHREELANDAULEVELSHAMILTONIAN_H
#define PARTICLEONTORUSTHREELANDAULEVELSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/BosonOnTorusWithSU3Spin.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnTorusThreeLandauLevelsHamiltonian : public ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian
{

 protected:

 //Parameters of the unit cell
 double Ratio;
 double InvRatio;
 double Lx;
 double Ly;
 double Gx;
 double Gy;

 //Number of fluxes
 int NbrLzValue;

 //Cyclotron energy
 double CyclotronEnergy;

 //Coulomb or delta interaction
 // true for delta, otherwise Coulomb
 bool HaveDelta;

 //Finite q=0 component in case of delta interaction
 double FiniteQZeroComponent;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = maximum Lz value reached by a particle in the state
  // ratio = ratio between the width in the x direction and the width in the y direction
  // haveDelta = true for delta interaction, otherwise Coulomb
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusThreeLandauLevelsHamiltonian(ParticleOnSphereWithSU3Spin* particles, int nbrParticles, int maxMomentum, double cyclotronenergy, double ratio, bool haveDelta,
						   AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusThreeLandauLevelsHamiltonian();

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for identical spins
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int l1, int l2, int l3, int l4);

  // Returns the Potential times the form factor
  // Q2 = momentum squared (Qx^2+Qy^2)
  // Qx = qx component
  // Qy = qy component 
  // l1,...,l4 = LL index
  Complex Potential(double Q2, double Qx, double Qy, int ll1, int ll2, int ll3, int ll4);

  // Returns the formfactor
  // Input : momentum + LL index
  Complex FormFactor(double Q2, double Qx, double Qy, int ll1, int ll2);
};

#endif
