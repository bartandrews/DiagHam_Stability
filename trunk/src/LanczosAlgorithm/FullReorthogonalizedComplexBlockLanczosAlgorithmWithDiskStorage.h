////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of full reorthogonalized complex block Lanczos algorithm     //
//                             using disk storage                             //
//                (with full re-orthogonalization at each step)               //
//                                                                            //
//                        last modification : 27/10/2014                      //
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


#ifndef FULLREORTHOGONALIZEDCOMPLEXBLOCKLANCZOSALGORITHMWITHDISKSTORAGE_H
#define FULLREORTHOGONALIZEDCOMPLEXBLOCKLANCZOSALGORITHMWITHDISKSTORAGE_H


#include "config.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexBlockLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/BandDiagonalHermitianMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"


class FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage : public FullReorthogonalizedComplexBlockLanczosAlgorithm
{

 protected:

  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

 public:

  // basic constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
  // blockSize = size of the block used for the block Lanczos algorithm
  // maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  // lapackFlag = rely on LAPACK library to diagonalize the block matrix
  FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize = 2, int maxIter = 1000,
								  bool strongConvergence = false, bool lapackFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage(const FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage& algorithm);

  // destructor
  //
  ~FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage();

  // get last produced vector
  //
  // return value = reference on last produced vector
  virtual Vector& GetGroundState();

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  virtual Vector* GetEigenstates(int nbrEigenstates);

  // resume Lanczos algorithm from disk datas in current directory
  //
  virtual void ResumeLanczosAlgorithm();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  virtual void RunLanczosAlgorithm (int nbrIter);
  
 protected:

  // reorthogonalize a set of vectors using Gram-Schmidt algorithm
  //
  // vectors = array of vectors to reorthogonalize
  // nbrVectors = number of vectors to reorthogonalize
  // matrix = matrix where transformation matrix has to be stored
  // rowShift = shift to apply to matrix row index to reach the upper leftmost element
  // columnShift = shift to apply to matrix column index to reach the upper leftmost element
  virtual void ReorthogonalizeVectors (ComplexVector* vectors, int nbrVectors, HermitianMatrix& matrix,
				       int rowShift, int columnShift);

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  virtual bool WriteState();

  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  virtual bool ReadState();

  // read a group of Lanczos vectors from disk
  // 
  // vectorArray = array where the vectors will be stored
  // vectorAbsoluteIndex = absolute index of the first vector that has to be read from disk
  // vectorRelativeIndex = index of the first vector that has to be read from disk within vectorArray
  // nbrVectors = number of vectors to read from disk
  // return value = true if no error occured  
  virtual bool ReadLanczosVectors(ComplexVector* vectorArray, int vectorAbsoluteIndex, int vectorRelativeIndex, int nbrVectors);

  // write a group of Lanczos vectors to disk
  // 
  // vectorArray = array where the vectors are stored
  // vectorAbsoluteIndex = absolute index of the first vector that has to be written to disk
  // vectorRelativeIndex = index of the first vector that has to be written to disk within vectorArray
  // nbrVectors = number of vectors to written to disk
  // return value = true if no error occured  
  virtual bool WriteLanczosVectors(ComplexVector* vectorArray, int vectorAbsoluteIndex, int vectorRelativeIndex, int nbrVectors);

};

#endif
