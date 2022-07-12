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


#ifndef MAXIMALLYCONDENSEDSTATEONLATTICE_H
#define MAXIMALLYCONDENSEDSTATEONLATTICE_H


#include "config.h"

#include "HilbertSpace/ParticleOnLattice.h"
#include "MathTools/Complex.h"
#include "MathTools/NSphereParameters.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"

#include <iostream>


using std::ostream;

class AbstractRandomNumberGenerator;

class MaximallyCondensedStateOnLattice
{
 protected:
  // number of sites
  int NbrVectors;
  
  // basis vectors
  ComplexVector* Vectors;

  // lattice dimensions
  int Lx;
  int Ly;
  int NbrSubLattices;

  // last value of largest density matrix eigenvalue
  double LastMaximumEV;

  // Hilbert space for vectors
  ParticleOnLattice *Space;

  // parametrization of variational parameters
  NSphereParameters SphereParametrization;
  
  // variational parameters
  RealVector VariationalParameters;

  // coordinates
  ComplexVector ResultingParameters;

  // number of eigenvalues to sum up
  int NbrEigenvalues;
  
  // number of evaluations in optimization
  int NbrEvaluations;

  // dimension of density matrices
  int DensityMatrixDimension;

  // storage for density matrices
  HermitianMatrix *DiagonalDensityMatrices;
  int NbrOffDiagonal;
  ComplexMatrix *OffDiagonalDensityMatrices;

  // current density matrix
  ComplexMatrix CurrentDensityMatrix;
  HermitianMatrix CurrentHermitianMatrix;
  
  // matrices as temporary space for calculations
  RealDiagonalMatrix M;
  ComplexMatrix Q;  

  // random number generator
  AbstractRandomNumberGenerator *RandomNumbers;

  // flag for external random number generator
  bool ExternalGenerator;
  
  // architecture for parallelization
  AbstractArchitecture *Architecture;
  
 public:

  // constructor for contact interactions on a square lattice
  //
  // nbrStates = number of quantum states
  // states = state vectors to superpose
  // space = Hilbert-space of states
  // lx = Lx dimension
  // ly = Ly dimension
  // sublattices = number of sublattices
  // randomGenerator = external random number generator
  MaximallyCondensedStateOnLattice(AbstractArchitecture *architecture, int nbrStates, ComplexVector *states, ParticleOnLattice* space, int lx, int ly, int sublattices, AbstractRandomNumberGenerator *randomGenerator=NULL);

  // destructor
  //
  ~MaximallyCondensedStateOnLattice();

  // get the parameters of the Many-Body state that was last calculated
  // return = state
  ComplexVector & GetVariationalParameters();

  // get the wavefunction corresponding to the current parameters
  // return = complex vector of local amplitudes and phases
  ComplexVector GetWaveFunction();

  // get last density matrix eigenvalue
  double GetDensityMatrixEigenvalue(){return this->LastMaximumEV;}

  // randomize trial parameters
  void RandomizeVariationalParameters();
  
  // set trial parameters
  void SetVariationalParameters(RealVector &variationalParameters);
  
  // optimize wavefunction starting from present settings of VariationalParameters
  // nbrEigenvals = number of eigenvalues that should be summed up for the optimization
  // tolerance = final tolerance on the variational parameters
  // maxIter = maximal number of function evaluations
  //
  double Optimize(int nbrEigenvals, double tolerance, int maxIter);

  double SimplexOptimize(double targetSize, int maxIter, double initialStep=1.0);

 private:

  // target function for optimizer routine:
  double EvaluateCondensateFraction(int nbrParameters, double *x);
  
  // evaluate all interaction factors
  //   
  void EvaluateDensityMatrices();
  
};

#endif // MAXIMALLYCONDENSEDSTATEONLATTICE_H
