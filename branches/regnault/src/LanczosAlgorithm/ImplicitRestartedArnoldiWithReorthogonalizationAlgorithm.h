////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of implicit restarted Arnoldi algorithm             //
//                          (with re-orthogonalization)                       //
//                                                                            //
//                        last modification : 06/01/2003                      //
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


#ifndef IMPLICITRESTARTEDARNOLDIWITHREORTHOGONALIZATIONALGORITHM_H
#define IMPLICITRESTARTEDARNOLDIWITHREORTHOGONALIZATIONALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/GarbageFlag.h"


class ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm : public AbstractLanczosAlgorithm
{

 protected:

  // array to store vectors of the invariant subspace 
  RealMatrix LanczosVectors;

  // number of needed eigenvectors
  int NbrEigenvectors;

  // number of Lanczos iteration per Arnoldi iteration
  int NbrIteration;

  // internal index corresponding to the number of implicit restarted Arnoldi steps realized
  int Index;

  // norm of the last vector produced during the standard Arnoldi step
  double LastVectorNorm;

  // vector to store store ground state
  RealVector GroundState;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // number of unwanted eigenvalues at each restarting
  int NbrUnwantedEigenvalues;
  // temporary array where unwanted eigenvalues are stored
  double* UnwantedEigenvalues;

  // error on the last wanted eigenvalue compare to its current associated Ritz value 
  double LastErrorRitzValue;

  double* ConvergedValues;
  int NbrConvergedValue;
  RealMatrix ConvergedEigenvectors;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvectors = number of needed eigenvectors (and eigenvalues)
  // nbrIteration = number of Lanczos iteration per Arnoldi iteration (equal 2 * nbrEigenvectors)
  ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm(AbstractArchitecture* architecture, int nbrEigenvectors, int nbrIteration);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm(const ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm& algorithm);

  // destructor
  //
  ~ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm();

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

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  Vector* GetEigenstates(int nbrEigenstates);

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();

 private:

  // standard Arnoldi step with full-reothogonalization
  //
  void StandardArnoldiStep();


};

#endif
