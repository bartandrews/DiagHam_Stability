////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of basic block Lanczos algorithm                   //
//                      (without re-orthogonalization )                       //
//                                                                            //
//                        last modification : 05/06/2008                      //
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


#ifndef BASICBLOCKLANCZOSALGORITHM_H
#define BASICBLOCKLANCZOSALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealBandDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class BasicBlockLanczosAlgorithm : public AbstractLanczosAlgorithm
{

 protected:

  // array where Lanczos vectors are stored
  RealVector* LanczosVectors;

  // vector array that contains the initial states
  RealVector* InitialStates;

  // maximum number of block Lanczos iteration
  int MaximumNumberIteration;

  // internal index corresponding to the number of block Lanczos steps that have been realized
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // vector to store store ground state
  RealVector GroundState;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;
  // value of the wanted eigenvalue at previous Lanczos iteration
  double* PreviousWantedEigenvalues;
  // flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  bool StrongConvergenceFlag;

  // size of the block used for the block Lanczos algorithm
  int BlockSize;

  // matrix where the blocks are stored
  RealBandDiagonalSymmetricMatrix ReducedMatrix;

  // temporary matrix used to duplicated ReducedMatrix before diagonalize it
  RealBandDiagonalSymmetricMatrix TemporaryReducedMatrix;
//  RealSymmetricMatrix TemporaryReducedMatrix;

  // array used to store temporary scalar products
  double* TemporaryCoefficients;

  // rely on LAPACK library to diagonalize the block matrix
  bool LapackFlag;

  // flag to indicate the use of disk storage to increase speed of ground state calculation
  bool DiskFlag;
  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
  // blockSize = size of the block used for the block Lanczos algorithm
  // maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  // lapackFlag = rely on LAPACK library to diagonalize the block matrix
  BasicBlockLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize = 2, int maxIter = 1000,
			     bool diskFlag  = false, bool resumeDiskFlag = false, bool strongConvergence = false, bool lapackFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicBlockLanczosAlgorithm(const BasicBlockLanczosAlgorithm& algorithm);

  // destructor
  //
  ~BasicBlockLanczosAlgorithm();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // initialize Lanczos algorithm with a set of given vectors
  //
  // vectors = array of vectors used as first step vectors
  // nbrVectors = number of vectors in the array
  void InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors);

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


 protected:
  
  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  bool ReadState();

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  bool WriteState();

  // diagonalize tridiagonalized matrix and find ground state energy
  //
  void Diagonalize ();

  // reorthogonalize a set of vectors using Gram-Schmidt algorithm
  //
  // vectors = array of vectors to reorthogonalize
  // nbrVectors = number of vectors to reorthogonalize
  // matrix = matrix where transformation matrix has to be stored
  // rowShift = shift to apply to matrix row index to reach the upper leftmost element
  // columnShift = shift to apply to matrix column index to reach the upper leftmost element
  void ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealBandDiagonalSymmetricMatrix& matrix,
                               int rowShift, int columnShift);
//  void ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealSymmetricMatrix& matrix,
//			       int rowShift, int columnShift);

};

#endif
