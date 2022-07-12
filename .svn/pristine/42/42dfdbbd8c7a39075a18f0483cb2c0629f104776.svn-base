////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of full reorthogonalized complex block Lanczos algorithm     //
//                (with full re-orthogonalization at each step)               //
//                                                                            //
//                        last modification : 05/01/2012                      //
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


#ifndef FULLREORTHOGONALIZEDCOMPLEXBLOCKLANCZOSALGORITHM_H
#define FULLREORTHOGONALIZEDCOMPLEXBLOCKLANCZOSALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/BandDiagonalHermitianMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"


class FullReorthogonalizedComplexBlockLanczosAlgorithm : public AbstractLanczosAlgorithm
{

 protected:

  // array where Lanczos vectors are stored
  ComplexVector* LanczosVectors;

  // maximum number of block Lanczos iteration
  int MaximumNumberIteration;

  // internal index corresponding to the number of block Lanczos steps that have been realized
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // vector to store store ground state
  ComplexVector GroundState;

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
//  BandDiagonalHermitianMatrix ReducedMatrix;
  HermitianMatrix ReducedMatrix;

  // temporary matrix used to duplicated ReducedMatrix before diagonalize it
//  BandDiagonalHermitianMatrix TemporaryReducedMatrix;
  HermitianMatrix TemporaryReducedMatrix;

  // array used to store temporary scalar products
  Complex* TemporaryCoefficients;

  // rely on LAPACK library to diagonalize the block matrix
  bool LapackFlag;

 public:

  // default constructor
  //
  FullReorthogonalizedComplexBlockLanczosAlgorithm();

  // basic constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
  // blockSize = size of the block used for the block Lanczos algorithm
  // maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  // lapackFlag = rely on LAPACK library to diagonalize the block matrix
  FullReorthogonalizedComplexBlockLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize = 2, int maxIter = 1000,
						   bool strongConvergence = false, bool lapackFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  FullReorthogonalizedComplexBlockLanczosAlgorithm(const FullReorthogonalizedComplexBlockLanczosAlgorithm& algorithm);

  // destructor
  //
  ~FullReorthogonalizedComplexBlockLanczosAlgorithm();

  // initialize Lanczos algorithm with a random vector
  //
  virtual void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  virtual void InitializeLanczosAlgorithm(const Vector& vector);

  // initialize Lanczos algorithm with a set of given vectors
  //
  // vectors = array of vectors used as first step vectors
  // nbrVectors = number of vectors in the array
  virtual void InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors);

  // get last produced vector
  //
  // return value = reference on last produced vector
  virtual Vector& GetGroundState();

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  virtual Vector* GetEigenstates(int nbrEigenstates);

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  virtual void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  virtual bool TestConvergence ();


 protected:
  
  // diagonalize tridiagonalized matrix and find ground state energy
  //
  virtual void Diagonalize ();


  // reorthogonalize a set of vectors using Gram-Schmidt algorithm
  //
  // vectors = array of vectors to reorthogonalize
  // nbrVectors = number of vectors to reorthogonalize
  // matrix = matrix where transformation matrix has to be stored
  // rowShift = shift to apply to matrix row index to reach the upper leftmost element
  // columnShift = shift to apply to matrix column index to reach the upper leftmost element
  //  void ReorthogonalizeVectors (ComplexVector* vectors, int nbrVectors, BandDiagonalHermitianMatrix& matrix,
  //                               int rowShift, int columnShift);
  virtual void ReorthogonalizeVectors (ComplexVector* vectors, int nbrVectors, HermitianMatrix& matrix,
			       int rowShift, int columnShift);


  virtual void TestOrthogonality (ComplexVector* vectors, int nbrVectors, ComplexVector* otherVectors = 0, int nbrOtherVectors = 0);

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
