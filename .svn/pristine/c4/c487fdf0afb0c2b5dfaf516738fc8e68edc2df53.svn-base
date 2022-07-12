////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of full reorthogonalized Lanczos algorithm             //
//                (with full re-orthogonalization at each step)               //
//                 and storing each iteration information on disk             //
//                                                                            //
//                        last modification : 18/03/2003                      //
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


#ifndef FULLREORTHOGONALIZEDLANCZOSALGORITHMWITHDISKSTORAGE_H
#define FULLREORTHOGONALIZEDLANCZOSALGORITHMWITHDISKSTORAGE_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class FullReorthogonalizedLanczosAlgorithmWithDiskStorage : public AbstractLanczosAlgorithm
{

 protected:

  RealVector* LanczosVectors;
  int MaximumNumberIteration;
  int Index;
  GarbageFlag Flag;
  RealVector GroundState;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;
  // value of the wanted eigenvalue at previous Lanczos iteration
  double* PreviousWantedEigenvalues;
  // flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  bool StrongConvergenceFlag;

  // number of vector that can be allocated at the same time in local memory (null if all vectors have to be stored in local memory)
  int MaxNbrVectors;

  // internally stored dimension, to be used if re-read from disk for vector calculation only
  int VectorDimension;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxNbrVectors = number of vector that can be allocated at the same time in local memory (null if all vectors have to be stored in local memory)
  // maxIter = an approximation of maximal number of iteration
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  FullReorthogonalizedLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int maxNbrVectors = 0, int maxIter = 100,
						      bool strongConvergence = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  FullReorthogonalizedLanczosAlgorithmWithDiskStorage(const FullReorthogonalizedLanczosAlgorithmWithDiskStorage& algorithm);

  // destructor
  //
  ~FullReorthogonalizedLanczosAlgorithmWithDiskStorage();

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
