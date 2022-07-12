////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2016 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quantum Hall hamiltonian associated               //
//       with the Kapit Mueller Hamiltonian for two coupled layers            //
//                with contact interactions on a lattice                      //
//                                                                            //
//                      last modification : 13/10/2016                        //
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


#ifndef PARTICLEONLATTICEPROJECTEDKAPITMUELLERMULTILAYERHAMILTONIAN_H
#define PARTICLEONLATTICEPROJECTEDKAPITMUELLERMULTILAYERHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>
#include <cmath>


using std::ostream;


class MathematicaOutput;
class KMBranchCut;


class ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian : public AbstractQHEOnLatticeHamiltonian
{
 protected:

  // number of (identical) layers
  int NbrLayers;

  // number of single particle states on which to project
  int NbrProjectorStates;

  // strength of intra-layer on-site delta-interaction
  double ContactInteractionU;

  // strength of inter-layer on-site delta-interaction
  double ContactInteractionW;

  // strength of delta-potential at the origin
  double DeltaPotential;

  // strength of random potential at individual sites
  double RandomPotential;

  // maximum range of single-particle hopping
  double Range;

  // parameters to pass on to EvaluateInteractionFactors:
  int NbrBranchCuts;
  double *BranchCoordinates;
  int *BranchShift;

  // One body eigenstate basis associated to each point of the band structure, the array index corresponds to the linearized momentum
  ComplexMatrix OneBodyBasis;

  // energy spectrum of the single-particle problem; index is the state index (no momentum implemented)
  double* EnergyLevels;

  // flag indicating whether to apply flat-band projection
  bool FlatBand;

  // flag for reversed hopping
  bool ReverseHopping;

 public:


  // constructor for contact interactions on a square lattice
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // numStates = number of single particle states to consider (projection to lowest numState SP states)
  // flatBand = apply the flat-band projection for all kept single-particle states
  // lx = length of simulation cell in x-direction
  // ly = length of simulation cell in y-direction
  // nbrLayers = number of layers to simulate
  // nbrFluxQuanta = number of flux quanta piercing the simulation cell
  // contactInteractionU = strength of intra-layer on-site delta interaction
  // contactInteractionW = strength of inter-layer on-site delta interaction
  // reverseHopping = flag to indicate if sign of hopping terms should be reversed
  // deltaPotential = strength of a delta potential at site (0,0)
  // randomPotential = magnitude of random potential to add to all sites
  // range = maximum length cut-off for single-particle hoppings
  // nbrBranchCuts = number of branch cuts
  // branchCoordinates coordinates of branch cuts in successive blocks of order (x_A, y_A, x_B, y_B)_i, i=1..nbrBranchCuts
  // branchShift shift of the layer index at the i-th branch cut
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // hermitianFlag = flag indicating whether to use hermitian symmetry
  ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian(ParticleOnLattice* particles, int nbrParticles, int numStates, bool flatBand, int lx, int ly, int nbrLayers, int nbrFluxQuanta, double contactInteractionU, double contactInteractionW, bool reverseHopping, double deltaPotential, double randomPotential, double range, int nbrBranchCuts, double *branchCoordinates, int *branchShift, AbstractArchitecture* architecture, int nbrBody = 2, unsigned long memory = 0, char* precalculationFileName = 0, bool hermitianFlag = false);

  // destructor
  //
  ~ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  
  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian& H);

 private:
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();
 
};


#endif // PARTICLEONLATTICEPROJECTEDKAPITMUELLERMULTILAYERHAMILTONIAN_H
