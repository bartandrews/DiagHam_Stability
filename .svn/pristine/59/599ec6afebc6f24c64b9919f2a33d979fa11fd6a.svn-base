////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2006 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of basic Lanczos algorithm with real vectors           //
//            and ground state evaluation and disk storage capabilities       //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 10/01/2006                      //
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


#ifndef BASICLANCZOSALGORITHMWITHGROUNDSTATEDISKSTORAGE_H
#define BASICLANCZOSALGORITHMWITHGROUNDSTATEDISKSTORAGE_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"


class BasicLanczosAlgorithmWithGroundStateDiskStorage : public AbstractLanczosAlgorithm
{

 protected:

  RealVector V1;
  RealVector V2;
  RealVector V3;

  // vector that contains the initial state
  RealVector InitialState;

  //vector that contains the ground state (if GroundStateFlag is set to true)
  RealVector GroundState;

  // flag to indicate if ground state has already been compute
  bool GroundStateFlag;

  // index to the current Lanczos iteration step  
  int Index;
  // flag to indicate if Index is related to the ground state calculation (1) or to the lowest eigenvalue calculation (0)
  int GroundStateEvaluationFlag;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;

  // maximum number of iterations when evaluating the ground state eigenvector 
  int NbrIterationsGroundState;

  // internally stored dimension, to be used if re-read from disk for vector calculation only
  int VectorDimension;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrIter = maximum number of iterations when evaluating the ground state eigenvector (0 if all iterations needed for convergence have to be done)
  // maxIter = an approximation of maximal number of iteration
  BasicLanczosAlgorithmWithGroundStateDiskStorage(AbstractArchitecture* architecture, int nbrIter = 0, int maxIter = 0);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicLanczosAlgorithmWithGroundStateDiskStorage(const BasicLanczosAlgorithmWithGroundStateDiskStorage& algorithm);

  // destructor
  //
  ~BasicLanczosAlgorithmWithGroundStateDiskStorage();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // resume Lanczos algorithm from disk datas in current directory
  //
  void ResumeLanczosAlgorithm();
  
  // get the n first eigenstates (limited to the ground state fro this class, return NULL if nbrEigenstates > 1)
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  Vector* GetEigenstates(int nbrEigenstates);

  // get ground state (by re-running Lanczos algorithm)
  //
  // return value = reference on ground state
  Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();
  
 private:

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  bool WriteState();

  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  bool ReadState();

};

#endif
