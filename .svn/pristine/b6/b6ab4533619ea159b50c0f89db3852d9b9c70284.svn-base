////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quantum Hall hamiltonian associated               //
//       with the Kapit Mueller Hamiltonian for two coupled layers            //
//                with contact interactions on a lattice                      //
//                                                                            //
//                      last modification : 13/04/2016                        //
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


#ifndef PARTICLEONLATTICEKAPITMUELLERBILAYERHAMILTONIAN_H
#define PARTICLEONLATTICEKAPITMUELLERBILAYERHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"

#include <iostream>
#include <cmath>


using std::ostream;


class MathematicaOutput;
class KMBranchCut;


class ParticleOnLatticeKapitMuellerMultiLayerHamiltonian : public AbstractQHEOnLatticeHamiltonian
{
 protected:

  // number of (identical) layers
  int NbrLayers;

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

  // flag for reversed hopping
  bool ReverseHopping;

  // number of branch cuts
  int NbrBranchCuts;

  // array with information about branch cuts
  KMBranchCut *BranchCuts;

  // magnitude of the shift in branch cuts
  int *BranchShift;


 public:


  // constructor for contact interactions on a square lattice
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lx = length of simulation cell in x-direction
  // ly = length of simulation cell in y-direction
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
  ParticleOnLatticeKapitMuellerMultiLayerHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, double contactInteractionW, bool reverseHopping, double deltaPotential, double randomPotential, double range, int nbrBranchCuts, double *branchCoordinates, int *branchShift, AbstractArchitecture* architecture, int nbrBody = 2, unsigned long memory = 0, char* precalculationFileName = 0, bool hermitianFlag = false);

  // destructor
  //
  ~ParticleOnLatticeKapitMuellerMultiLayerHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  
  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnLatticeKapitMuellerMultiLayerHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeKapitMuellerMultiLayerHamiltonian& H);

 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  /// evaluate branch crossings for a link from site i to j
  /// @param xi, yi, xj, yj coordinates of the initial (i) and final site (j)
  /// @return overall shift for crossing any of the known branch cuts
  int EvaluateBranchCrossings(int xi, int yi, int xj, int yj);

  
  // calculate hopping amplitude and phase for the Kapit-Mueller single-particle Hamiltonian
  Complex KapitMuellerHopping(int xj, int yj, int xi, int yi);

  Complex SumImagesForHoppings(int xj, int yj, int xi, int yi, int limit);

  // calculate distance on lattice
  double GetDistance(int dx, int dy) {return sqrt(dx*dx+dy*dy);}

};


inline Complex ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::KapitMuellerHopping(int xj, int yj, int xi, int yi)
{
  int x = xj - xi;
  int y = yj - yi;
  int sgn = 1-2*((x + y + x*y) & 0x1);

  double amplitude = sgn*std::exp(-0.5*M_PI*(1.0-this->FluxDensity)*((x*x + y*y)));
  if (this->Particles->GetLandauGaugeAxis()=='y')
    return Polar(amplitude, -M_PI * this->FluxDensity * (yj-yi) * (xj+xi) );
  else if (this->Particles->GetLandauGaugeAxis()=='x')
    return Polar(amplitude, -M_PI * this->FluxDensity *  (yj+yi) * (xj-xi) );
  else
    {
      std::cerr << "Unknown Landau quantization axis"<<endl;
      exit(1);
    }
}


inline Complex ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::SumImagesForHoppings(int xj, int yj, int xi, int yi, int images)
{
  Complex sum=0.0;
  double GaugeX, GaugeY;
  this->Particles->GetSolenoidFluxes(GaugeX, GaugeY);
  for (int dX = -images; dX <= images; ++dX)
    for (int dY = -images; dY <= images; ++dY)
      {
	sum+=Polar(-2.0*M_PI*this->FluxDensity*this->Lx*dX * yj + dX*GaugeX + dY*GaugeY) * this->KapitMuellerHopping(xj+dX*this->Lx, yj+dY*this->Ly, xi, yi);
      }
  std::cout << "JL( ("<<xj<<", "<<yj<<")<-("<<xi<<", "<<yi<<")="<<sum<<std::endl;
  return sum;
}

#endif // PARTICLEONLATTICEKAPITMUELLERBILAYERHAMILTONIAN_H
