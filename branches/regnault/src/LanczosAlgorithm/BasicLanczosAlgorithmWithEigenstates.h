////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of basic Lanczos algorithm                     //
//                          with access to eigenstates                        //
//                                                                            //
//                        last modification : 17/07/2001                      //
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


#ifndef BASICLANCZOSALGORITHMWITHEIGENSTATES_H
#define BASICLANCZOSALGORITHMWITHEIGENSTATES_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class BasicLanczosAlgorithmWithEigenstates : public AbstractLanczosAlgorithm
{

 protected:

  RealVector* LanczosVectors;
  int MaximumNumberIteration;
  int Index;
  GarbageFlag Flag;
  RealVector GroundState;

 public:

  // default constructor
  //
  // maxIter = an approximation of maximal number of iteration
  BasicLanczosAlgorithmWithEigenstates(AbstractArchitecture* architecture, int maxIter = 100);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicLanczosAlgorithmWithEigenstates(const BasicLanczosAlgorithmWithEigenstates& algorithm);

  // destructor
  //
  ~BasicLanczosAlgorithmWithEigenstates();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // get ground state
  //
  // return value = reference on ground state
  Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
};

#endif
