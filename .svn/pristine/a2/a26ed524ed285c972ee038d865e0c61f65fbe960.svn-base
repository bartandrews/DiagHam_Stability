////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of power methof algorithm                      //
//                                                                            //
//                        last modification : 24/01/2013                      //
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


#ifndef POWERMETHODALGORITHM_H
#define POWERMETHODALGORITHM_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/RealUpperHessenbergMatrix.h"


class PowerMethodAlgorithm : public AbstractLanczosAlgorithm
{

 protected:


  // array where the vectors of the Krylov subspace
  RealVector PreviousVector;

  // array where the vectors of the Krylov subspace
  RealVector CurrentVector;

  // maximum  number of iterations
  int MaximumNumberIteration;
  // current iteration index
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // value of the last wanted eigenvalue at previous Arnoldi iteration
  double PreviousEigenvalue;
  // value of the wanted eigenvalue at previous Arnoldi iteration
  double CurrentEigenvalue;

 public:

  // default constructor
  //
  PowerMethodAlgorithm();

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  PowerMethodAlgorithm(AbstractArchitecture* architecture, int maxIter = 100);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  PowerMethodAlgorithm(const PowerMethodAlgorithm& algorithm);

  // destructor
  //
  ~PowerMethodAlgorithm();

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

  // get the n first eigenvalues
  //
  // eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
  // nbrEigenstates = number of needed eigenvalues
  virtual void GetEigenvalues (double*& eigenvalues, int nbrEigenvalues);

  // get the n first eigenvalues
  //
  // eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
  // nbrEigenstates = number of needed eigenvalues
  virtual void GetEigenvalues (Complex*& eigenvalues, int nbrEigenvalues);

  // run current Arnoldi algorithm (continue from previous results if Arnoldi algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();

};

#endif
