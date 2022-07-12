////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                                                                            //
//                        last modification : 23/02/2011                      //
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


#ifndef PARTICLEONLATTICECHERNINSULATORSINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICECHERNINSULATORSINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;
class Abstract2DTightBindingModel;
class Polynomial;


class ParticleOnLatticeChernInsulatorSingleBandHamiltonian : public ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian
{

 protected:
  
  // pointer to the tight binding model
  Abstract2DTightBindingModel* TightBindingModel;

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;

  // index of the band that has to be partially filled
  int BandIndex;

  // the interpolation parameter lambda between FQH and FCI, as in H = lambda * H_FQH + (1 - lambda) * H_FCI
  double InterpolationToFQH;
  // dimension of the internal (color) space. this should be equal to the Chern number of the occupid band
  int NbrColor;
  // angle between the two fundamental cycles of the torus in Radians
  double TwistAngle;
  // aspect ratio of torus, Lx / Ly
  double AspectRatio;
  // color-entangled LLL boundary condition twisting angle along x
  double LLLGammaX;
  // color-entangled LLL boundary condition twisting angle along y
  double LLLGammaY;
  // array of the pseudo-potentials
  double* Pseudopotentials;
  // array of the number of pseudo-potentials
  int NbrPseudopotentials;
  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;

  
 public:

  // default constructor
  //
  ParticleOnLatticeChernInsulatorSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tightBindingModel = pointer to the tight binding model
  // nbrColor = dimension of the internal (color) space. this should be equal to the Chern number of the occupid band
  // twistAngle = angle between the two fundamental cycles of the torus in Radians
  // aspectRatio = aspect ratio of torus, Lx / Ly
  // lLLGammaX = color-entangled LLL boundary condition twisting angle along x
  // lLLGammaY = color-entangled LLL boundary condition twisting angle along y
  // nbrPseudopotentials = array of the number of pseudo-potentials
  // pseudoPotentials = array of the pseudo-potentials
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeChernInsulatorSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel,
          int nbrColor, double twistAngle, double aspectRatio, double lLLGammaX, double lLLGammaY, int nbrPseudopotentials, double* pseudoPotentials,
          AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeChernInsulatorSingleBandHamiltonian();
  

 protected:
 
  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices
  virtual void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // perform gauge transform to the LLL gauge
  //
  void TransformToLLLGauge(ComplexMatrix& gauge, double gammaX, double gammaY);

  // evaluate all FQH interaction factors
  //
  virtual void EvaluateFQHInteractionFactors();

  // evaluate matrix element of the FQH pseudopotential Hamiltonian, namely the numerical coefficient in front of the a+_k1 a+_k2 a_k3 a_k4
  // kx1 = first kx index
  // ky1 = first ky index
  // kx2 = second kx index
  // ky2 = second ky index
  // kx3 = third kx index
  // ky3 = third ky index
  // kx4 = fourth kx index
  // ky4 = fourth ky index
  // nbrPseudopotentials = number of pseudopotentials
  // pseudopotentials = pseudopotential coefficients
  // return value = numerical coefficient
  virtual Complex EvaluateFQHInteractionCoefficient(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // get fourier transform of FQH pseudopotential interaction times the exponential supression
  // Q2_half = one half of q² value
  double GetVofQ(double Q2_half);

};


#endif
