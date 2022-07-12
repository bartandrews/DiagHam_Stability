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


#ifndef EXPLICITMOOREREADONSPHEREWAVEFUNCTION_H
#define EXPLICITMOOREREADONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "MathTools/NumericalAnalysis/SymmetrizedComplexFunction.h"
#include "GeneralTools/GarbageFlag.h"

#include "ReadRezayiOnSphereBlock.h"


class ExplicitMooreReadOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
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

  // internal real coordinates
  RealVector RealCoordinates;

  SymmetrizedComplexFunction *Symmetrizer;

  ReadRezayiOnSphereBlock *Block;
  
  // garable flag associated to the Permutations array
  GarbageFlag Flag;

 public:
  
  // constructor
  //
  // nbrParticlesPerCluster = number of particles per cluster
  // nbrClusters = number of clusters
  // fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
  ExplicitMooreReadOnSphereWaveFunction(int nbrParticlesPerCluster, int nbrClusters, bool fermionicStatistics = true);

  // copy constructor
  //
  // function = reference on the wave function to copy
  ExplicitMooreReadOnSphereWaveFunction(const ExplicitMooreReadOnSphereWaveFunction& function);

  // destructor
  //
  ~ExplicitMooreReadOnSphereWaveFunction();

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


};

#endif
