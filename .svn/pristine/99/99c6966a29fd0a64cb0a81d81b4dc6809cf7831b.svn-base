////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of implicetly restarted Arnoldi algorithm              //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 06/02/2013                      //
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


#ifndef IMPLICITLYRESTARTEDARNOLDIALGORITHM_H
#define IMPLICITLYRESTARTEDARNOLDIALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/RealUpperHessenbergMatrix.h"


class ImplicitlyRestartedArnoldiAlgorithm : public BasicArnoldiAlgorithm
{

 protected:

  // maximum number of vectors that can be stored in memory before restarting the Arnoldi algorithm
  int MaxNbrVectors;
  // number of vectors that are kept when restarting the Arnoldi algorithm
  int NbrKeptVectors;

  // matrix where the Arnoldi vectors are stored (wrapper to ArnoldiVectors)
  RealMatrix ArnoldiVectorMatrix;

 public:

  // default constructor
  //
  ImplicitlyRestartedArnoldiAlgorithm();

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxNbrVectors = maximum number of vectors that can be stored in memory before restarting the Arnoldi algorithm
  // nbrKeptVectors = number of vectors that are kept when restarting the Arnoldi algorithm (can't be lower than nbrEigenvalue)
  // maxIter = an approximation of maximal number of iteration
  // highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  // leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  ImplicitlyRestartedArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxNbrVectors, int nbrKeptVectors,
				      int maxIter = 100, bool highEnergy = false, bool leftFlag = false, bool strongConvergence = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ImplicitlyRestartedArnoldiAlgorithm(const ImplicitlyRestartedArnoldiAlgorithm& algorithm);

  // destructor
  //
  ~ImplicitlyRestartedArnoldiAlgorithm();

  // initialize Arnoldi algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Arnoldi algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

 protected:

  // restart the Arnoldi algorithm if needed
  //
  virtual void RestartAlgorithm();


};

#endif
