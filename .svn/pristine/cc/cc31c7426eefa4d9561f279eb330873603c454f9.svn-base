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


#ifndef BASICLANCZOSALGORITHMWITHGROUNDSTATE_H
#define BASICLANCZOSALGORITHMWITHGROUNDSTATE_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"


class BasicLanczosAlgorithmWithGroundState : public AbstractLanczosAlgorithm
{

 protected:

  RealVector V1;
  RealVector V2;
  RealVector V3;

  // vector that contains the initial state
  RealVector InitialState;

  //vector that contains the ground state (if GroundStateFlag is set to true)
  RealVector GroundState;

  // flag to indicate if ground state has already been computed
  bool GroundStateFlag;

  // flag to indicate the use of disk storage to increase speed of ground state calculation
  bool DiskFlag;
  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

  int Index;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;

  // number of vectors to which eigenstates have to be othogonal
  int OrthogonalizationSetSize;
  // vectors to which eigenstates have to be othogonal
  RealVector* OrthogonalizationSet;
  // names of the files corresponding to vectors to which eigenstates have to be othogonal
  char** OrthogonalizationSetFileNames;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  BasicLanczosAlgorithmWithGroundState(AbstractArchitecture* architecture, int maxIter = 0, bool diskFlag = false, bool resumeDiskFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicLanczosAlgorithmWithGroundState(const BasicLanczosAlgorithmWithGroundState& algorithm);

  // destructor
  //
  ~BasicLanczosAlgorithmWithGroundState();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // force orthogonalization with respect to a set of vectors
  //
  // fileName = name of the file describing the set of vectors
  // return value = true if no error occured
  bool ForceOrthogonalization(char* fileName);

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
  
 protected:
  
  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  bool ReadState();

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  bool WriteState();

  // orthogonalize a vector with respect to a set of external vectors
  //
  // inputVector = reference on the vector whose component on the external set has to be removed 
  void ExternalOrthonogalization(RealVector& inputVector);

};

#endif
