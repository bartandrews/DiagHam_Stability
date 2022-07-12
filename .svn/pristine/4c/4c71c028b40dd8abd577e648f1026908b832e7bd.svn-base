////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of basic Lanczos algorithm with complex vectors         //
//                 and ground state evaluation, using both fast disk          //
//                     and a projector over a set of vectors                  //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 25/05/2015                      //
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


#ifndef COMPLEXBASICLANCZOSALGORITHMWITHGROUNDSTATEANDPROJECTORFASTDISK_H
#define COMPLEXBASICLANCZOSALGORITHMWITHGROUNDSTATEANDPROJECTORFASTDISK_H


#include "config.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/ComplexVector.h"


class ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk : public ComplexBasicLanczosAlgorithmWithGroundStateFastDisk
{

 protected:

  // dimension of the projector subspace
  int NbrProjectors;
  // initial dimension of the projector subspace (before it is increased by the automatic generation)
  int InitialNbrProjectors;

  // array that contains the vectors that spans the projector subspace
  ComplexVector* ProjectorVectors;
  // energy scale in front of the projector
  double ProjectorCoefficient;
  // true if the eigenstate indices have to be shifted
  bool IndexShiftFlag;

  //  true if the projector subspace has to be constructed automatically
  bool AutomaticProjectorConstructionFlag;
  // array that contains the eigenvalue of each state of the projector subspace (only useful when using the automatic projector subspace generation)
  double* ProjectorEigenvalues;
  // temporary matrix used for the automatic projector subspace generation
  RealTriDiagonalSymmetricMatrix FullDiagonalizedMatrix;

  // flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
  bool ReturnAtConvergenceFlag;
  
 public:

  // default constructor
  //
  // nbrProjectors = dimension of the projector subspace
  // projectorVectors = array that contains the vectors that spans the projector subspace
  // projectorCoefficient = energy scale in front of the projector
  // indexShiftFlag = true if the eigenstate indices have to be shifted
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  // returnAtConvergenceFlag = flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
  ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrProjectors, ComplexVector* projectorVectors, double projectorCoefficient, 
								  bool indexShiftFlag,
								  AbstractArchitecture* architecture, 
								  int maxIter = 0, bool diskFlag = false, bool resumeDiskFlag = false, bool returnAtConvergenceFlag = false);

  // constructor using automatic projector construction 
  //
  // nbrEigenvalues = number of eigenvalues/eigenstates to compute
  // projectorCoefficient = energy scale in front of the projector
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrEigenvalues, double projectorCoefficient,
								  AbstractArchitecture* architecture, 
								  int maxIter = 0, bool diskFlag = false, bool resumeDiskFlag = false, bool returnAtConvergenceFlag=false);

  // constructor using both automatic projector construction and an initial set of projectors
  //
  // nbrEigenvalues = number of eigenvalues/eigenstates to compute
  // nbrProjectors = dimension of the projector subspace
  // projectorVectors = array that contains the vectors that spans the projector subspace
  // projectorCoefficient = energy scale in front of the projector
  // indexShiftFlag = true if the eigenstate indices have to be shifted
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrEigenvalues, 
								  int nbrProjectors, ComplexVector* projectorVectors,
								  double projectorCoefficient, bool indexShiftFlag,
								  AbstractArchitecture* architecture, 
								  int maxIter = 0, bool diskFlag = false, bool resumeDiskFlag = false, bool returnAtConvergenceFlag=false);
  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(const ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk& algorithm);

  // destructor
  //
  ~ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk();

  // initialize Lanczos algorithm with a random vector
  //
  virtual void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  virtual void InitializeLanczosAlgorithm(const Vector& vector);

  // get the n first eigenstates (limited to the ground state fro this class, return NULL if nbrEigenstates > 1)
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  virtual Vector* GetEigenstates(int nbrEigenstates);

  // get ground state (by re-running Lanczos algorithm)
  //
  // return value = reference on ground state
  virtual Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  virtual void RunLanczosAlgorithm (int nbrIter);
  
  // optional shift of the eigenstate file name indices
  //
  // return value = index shift
  virtual int EigenstateIndexShift();

  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  virtual bool TestConvergence ();
  
  // get the n first eigenvalues
  //
  // eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
  // nbrEigenstates = number of needed eigenvalues
  virtual void GetEigenvalues (double*& eigenvalues, int nbrEigenvalues);

  // get current diagonalized matrix
  //
  // return value = reference on current diagonalized matrix
  virtual RealTriDiagonalSymmetricMatrix& GetDiagonalizedMatrix ();

 protected:
  
  // add the projector contribution to the hamiltonian-vector multiplication
  //
  // initialVector = reference on the initial vector
  // destinationVector = reference on the destination vector 
  virtual void AddProjectorContribution(ComplexVector& initialVector, ComplexVector& destinationVector);

  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  virtual bool ReadState();

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  virtual bool WriteState();

};

#endif
