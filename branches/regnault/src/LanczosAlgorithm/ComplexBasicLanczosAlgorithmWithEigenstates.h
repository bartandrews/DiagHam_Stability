////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of basic Lanczos algorithm with complex vectors          //
//                           and access to eigenstates                        //
//                                                                            //
//                        last modification : 26/03/2002                      //
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


#ifndef COMPLEXBASICLANCZOSALGORITHMWITHEIGENSTATES_H
#define COMPLEXBASICLANCZOSALGORITHMWITHEIGENSTATES_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"


class ComplexBasicLanczosAlgorithmWithEigenstates : public AbstractLanczosAlgorithm
{

 protected:

  ComplexVector* LanczosVectors;
  int MaximumNumberIteration;
  GarbageFlag Flag;
  ComplexVector GroundState;
  int Index;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  ComplexBasicLanczosAlgorithmWithEigenstates(AbstractArchitecture* architecture, int maxIter = 100);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ComplexBasicLanczosAlgorithmWithEigenstates(const ComplexBasicLanczosAlgorithmWithEigenstates& algorithm);

  // destructor
  //
  ~ComplexBasicLanczosAlgorithmWithEigenstates();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // get last produced vector
  //
  // return value = reference on last produced vector
  Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
};

#endif
