////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of basic Lanczos algorithm with complex vectors          //
//                     and a projector over a set of vectors                  //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 11/06/2015                      //
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


#ifndef COMPLEXBASICLANCZOSALGORITHMWITHPROJECTOR_H
#define COMPLEXBASICLANCZOSALGORITHMWITHPROJECTOR_H


#include "config.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/ComplexVector.h"


class ComplexBasicLanczosAlgorithmWithProjector : public ComplexBasicLanczosAlgorithm
{

 protected:

  // dimension of the projector subspace
  int NbrProjectors;
  // array that contains the vectors that spans the projector subspace
  ComplexVector* ProjectorVectors;
  // energy scale in front of the projector
  double ProjectorCoefficient;
  // true if the eigenstate indices have to be shifted
  bool IndexShiftFlag;

 public:

  // default constructor
  //
  // nbrProjectors = dimension of the projector subspace
  // projectorVectors = array that contains the vectors that spans the projector subspace
  // projectorCoefficient = energy scale in front of the projector
  // indexShiftFlag = true if the eigenstate indices have to be shifted
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxIter = an approximation of maximal number of iteration
  ComplexBasicLanczosAlgorithmWithProjector(int nbrProjectors, ComplexVector* projectorVectors, double projectorCoefficient, bool indexShiftFlag,
					    AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter = 0);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ComplexBasicLanczosAlgorithmWithProjector(const ComplexBasicLanczosAlgorithmWithProjector& algorithm);

  // destructor
  //
  ~ComplexBasicLanczosAlgorithmWithProjector();

  // initialize Lanczos algorithm with a random vector
  //
  virtual void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  virtual void InitializeLanczosAlgorithm(const Vector& vector);

  // get last produced vector
  //
  // return value = reference on last produced vector
  virtual Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  virtual void RunLanczosAlgorithm (int nbrIter);
  
 protected:
  
  // add the projector contribution to the hamiltonian-vector multiplication
  //
  // initialVector = reference on the initial vector
  // destinationVector = reference on the destination vector 
  virtual void AddProjectorContribution(ComplexVector& initialVector, ComplexVector& destinationVector);

};

#endif
