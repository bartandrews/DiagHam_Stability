////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                 for non symmetric matrices using disk storage              //
//                                                                            //
//                        last modification : 07/01/2013                      //
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


#ifndef BASICARNOLDIALGORITHMWITHDISKSTORAGE_H
#define BASICARNOLDIALGORITHMWITHDISKSTORAGE_H


#include "config.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/RealUpperHessenbergMatrix.h"


class BasicArnoldiAlgorithmWithDiskStorage : public BasicArnoldiAlgorithm
{

 protected:

  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

  // maximum number of vectors that can be stored
  int NbrTemporaryVectors;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxIter = an approximation of maximal number of iteration
  // highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  // leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  // nbrTemporaryVectors = number of temporary that can be stored in memory
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  BasicArnoldiAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter = 100,
				       bool highEnergy = false, bool leftFlag = false, bool resumeDiskFlag = false, int nbrTemporaryVectors = 1, bool strongConvergence = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicArnoldiAlgorithmWithDiskStorage(const BasicArnoldiAlgorithmWithDiskStorage& algorithm);

  // destructor
  //
  ~BasicArnoldiAlgorithmWithDiskStorage();

  // initialize Arnoldi algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Arnoldi algorithm with a given vector
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

  // run current Arnoldi algorithm (continue from previous results if Arnoldi algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
 protected:

  // read several temporary vectors stored o disk
  //
  // firstVector = index of the first vector to read
  // totalNbrVectors = total number of temporary vectors
  // return value = number of vectors that have been read
  int ReadTemporaryVectors(int firstVector, int totalNbrVectors);

  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  bool ReadState();

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  bool WriteState();


};

#endif
