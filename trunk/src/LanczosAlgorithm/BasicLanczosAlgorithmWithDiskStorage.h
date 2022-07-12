////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of basic Lanczos algorithm                     //
//                  (without re-orthogonalization at each step)               //
//                 and storing each iteration information on disk             //
//                                                                            //
//                        last modification : 16/03/2003                      //
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


#ifndef BASICLANCZOSALGORITHMWITHDISKSTORAGE_H
#define BASICLANCZOSALGORITHMWITHDISKSTORAGE_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"


class BasicLanczosAlgorithmWithDiskStorage : public AbstractLanczosAlgorithm
{

 protected:

  RealVector V1;
  RealVector V2;
  RealVector V3;

  int Index;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;

  // file name containing  Lanczos iteration informations
  char* LanczosFileName;

  // file name corresponding to the first Lanczos vector
  char* V1FileName;
  // file name corresponding to the second Lanczos vector
  char* V2FileName;
  // file name corresponding to the third Lanczos vector
  char* V3FileName;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxIter = an approximation of maximal number of iteration
  BasicLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter = 0);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  BasicLanczosAlgorithmWithDiskStorage(const BasicLanczosAlgorithmWithDiskStorage& algorithm);

  // destructor
  //
  ~BasicLanczosAlgorithmWithDiskStorage();

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
