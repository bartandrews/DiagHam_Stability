////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of basic block Arnoldi algorithm                   //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 05/12/2012                      //
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


#ifndef BASICBLOCKARNOLDIALGORITHM_H
#define BASICBLOCKARNOLDIALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class BasicBlockArnoldiAlgorithm : public AbstractLanczosAlgorithm
{

 protected:

  // array where the vectors of the Krylov subspace
  RealVector* ArnoldiVectors;

  // maximum  number of iterations
  int MaximumNumberIteration;
  // current iteration index
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // Hessenberg matrix of the Hamiltonian in the Krylov subspace
  RealMatrix ReducedMatrix;
  // temporary matrix used to duplicated ReducedMatrix before diagonalize it
  RealMatrix TemporaryReducedMatrix;

  ComplexDiagonalMatrix ComplexDiagonalizedMatrix;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Arnoldi iteration
  Complex PreviousLastWantedEigenvalue;
  // value of the wanted eigenvalue at previous Arnoldi iteration
  Complex* ComplexPreviousWantedEigenvalues;
  // ground state 
  ComplexVector GroundState;

  // flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  bool StrongConvergenceFlag;
  //true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  bool HighEnergyFlag;
  // flag to compute the left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  bool LeftFlag;

  // array used to store temporary scalar products
  double* TemporaryCoefficients;

  // size of the block used for the block Lanczos algorithm
  int BlockSize;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // blockSize = size of the block used for the block Lanczos algorithm
  // maxIter = an approximation of maximal number of iteration
  // highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  // leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  BasicBlockArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize = 2, int maxIter = 100,
			     bool highEnergy = false, bool leftFlag = false, bool strongConvergence = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicBlockArnoldiAlgorithm(const BasicBlockArnoldiAlgorithm& algorithm);

  // destructor
  //
  ~BasicBlockArnoldiAlgorithm();

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
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();

 protected:

  // diagonalize tridiagonalized matrix and find ground state energy
  //
  virtual void Diagonalize();

  // reorthogonalize a set of vectors using Gram-Schmidt algorithm
  //
  // vectors = array of vectors to reorthogonalize
  // nbrVectors = number of vectors to reorthogonalize
  // matrix = matrix where transformation matrix has to be stored
  // rowShift = shift to apply to matrix row index to reach the upper leftmost element
  // columnShift = shift to apply to matrix column index to reach the upper leftmost element
  virtual void ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealMatrix& matrix,
				       int rowShift, int columnShift);
};

#endif
