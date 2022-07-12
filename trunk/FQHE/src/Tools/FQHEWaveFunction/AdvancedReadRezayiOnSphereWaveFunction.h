////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Moore Read state wave function on sphere             //
//                                                                            //
//                        last modification : 19/09/2004                      //
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


#ifndef ADVANCEDREADREZAYIONSPHEREWAVEFUNCTION_H
#define ADVANCEDREADREZAYIONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "GeneralTools/GarbageFlag.h"

class AdvancedReadRezayiOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles
  int NbrParticles;

  // number of particle per cluster
  int ClusterSize;

  // number of clusters
  int NbrClusters;

  // flag indicating whether we have fermions
  bool FermionicStatistics;

  // internal storage of spinor coordinates
  Complex* SpinorUCoordinates;
  Complex* SpinorVCoordinates;

  // individual JastrowFactor terms (squared)
  Complex** JastrowFactorElements;
  Complex** JastrowFactorSquares;

  // Jastrow-Factors for Block k
  Complex* BlockJk;
  // array containing description of each permutation that appears in the calculation of the Moore-Read state
  unsigned** Permutations;
  // number of different permutations that appear in the calculation of the Moore-Read state
  unsigned NbrPermutations;
  // weight of permutations
  double* WeightOfPermutations;
  // garable flag associated to the Permutations array
  GarbageFlag Flag;
  

 public:
  
  // constructor
  //
  // nbrParticlesPerCluster = number of particles per cluster (=N/2)
  // nbrClusters = number of clusters
  // fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
  AdvancedReadRezayiOnSphereWaveFunction(int nbrParticlesPerCluster, int nbrClusters, bool fermionicStatistics = true);

  // copy constructor
  //
  // function = reference on the wave function to copy
  AdvancedReadRezayiOnSphereWaveFunction(const AdvancedReadRezayiOnSphereWaveFunction& function);

  // destructor
  //
  ~AdvancedReadRezayiOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

 private:

  // evaluate all permutations requested for the Moore-Read state evaluation
  //
  void EvaluatePermutations();

  // perform complex part of calculations
  // uses internal spinor coordinates as input
  //
  Complex ComplexEvaluations();

};

#endif
