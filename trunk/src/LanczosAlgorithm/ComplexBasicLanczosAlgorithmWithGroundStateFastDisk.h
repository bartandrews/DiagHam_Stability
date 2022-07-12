////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of basic Lanczos algorithm with real vectors           //
//                         and ground state evaluation                        //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 17/09/2002                      //
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


#ifndef COMPLEXBASICLANCZOSALGORITHMWITHGROUNDSTATEFASTDISK_H
#define COMPLEXBASICLANCZOSALGORITHMWITHGROUNDSTATEFASTDISK_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/ComplexVector.h"


class ComplexBasicLanczosAlgorithmWithGroundStateFastDisk : public AbstractLanczosAlgorithm
{

 protected:

  ComplexVector V1;
  ComplexVector V2;
  ComplexVector V3;

  // vector that contains the initial state
  ComplexVector InitialState;

  //vector that contains the ground state (if GroundStateFlag is set to true)
  ComplexVector GroundState;

  // flag to indicate if ground state has already been compute
  bool GroundStateFlag;

  // flag to indicate the use of disk storage to increase speed of ground state calculation
  bool DiskFlag;
  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

  // flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
  bool ReturnAtConvergenceFlag;

  int Index;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;

 public:

  // default constructor
  //
  ComplexBasicLanczosAlgorithmWithGroundStateFastDisk();

  // constructor
  //
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  // returnAtConvergenceFlag = flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
  ComplexBasicLanczosAlgorithmWithGroundStateFastDisk(AbstractArchitecture* architecture, int maxIter = 0, bool diskFlag = false, bool resumeDiskFlag = false, bool returnAtConvergenceFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ComplexBasicLanczosAlgorithmWithGroundStateFastDisk(const ComplexBasicLanczosAlgorithmWithGroundStateFastDisk& algorithm);

  // destructor
  //
  ~ComplexBasicLanczosAlgorithmWithGroundStateFastDisk();

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
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  virtual bool TestConvergence ();
  
 protected:
  
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
