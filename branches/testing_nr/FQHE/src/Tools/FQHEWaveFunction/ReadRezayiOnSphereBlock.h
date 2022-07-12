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


#ifndef READREZAYIONSPHEREBLOCK_H
#define READREZAYIONSPHEREBLOCK_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "GeneralTools/GarbageFlag.h"


class ReadRezayiOnSphereBlock: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // number of particle per cluster
  int ClusterSize;

  // number of clusters
  int NbrClusters;

  // internal storage of spinor coordinates
  Complex* SpinorUCoordinates;
  Complex* SpinorVCoordinates;

  // garable flag associated to the Permutations array
  GarbageFlag Flag;

 public:
  
  // constructor
  //
  // nbrParticlesPerCluster = number of particles per cluster (=N/2)
  // nbrCluster = number of clusters
  ReadRezayiOnSphereBlock(int nbrParticlesPerCluster, int nbrCluster);
    
  // copy constructor
  //
  // function = reference on the wave function to copy
  ReadRezayiOnSphereBlock(const ReadRezayiOnSphereBlock& function);

  // destructor
  //
  ~ReadRezayiOnSphereBlock();

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
  
};

#endif
