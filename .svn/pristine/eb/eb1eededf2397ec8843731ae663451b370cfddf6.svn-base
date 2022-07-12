////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of Moore Read state wave function on disk              //
//                                                                            //
//                        last modification : 05/11/2004                      //
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


#ifndef MOOREREADONDISKWAVEFUNCTION_H
#define MOOREREADONDISKWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "GeneralTools/GarbageFlag.h"


class MooreReadOnDiskWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // number of particle per cluster
  int ClusterSize;

  // number of clusters
  int NbrClusters;

  // array containing description of each permutation that appears in the calculation of the Moore-Read state
  unsigned long** Permutations;
  // number of permutations that appears in the calculation of the Moore-Read state
  unsigned long NbrPermutations;
  // garbage flag associated to the Permutations array
  GarbageFlag Flag;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // clusterSize = number of particle per cluster
  MooreReadOnDiskWaveFunction(int nbrParticles, int clusterSize);

  // copy constructor
  //
  // function = reference on the wave function to copy
  MooreReadOnDiskWaveFunction(const MooreReadOnDiskWaveFunction& function);

  // destructor
  //
  ~MooreReadOnDiskWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

 private:

  // evaluate all permutations requested for the Moore-Read state evaluation
  //
  void EvaluatePermutations();

};

#endif
