////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 17/11/2012                      //
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


#ifndef COMPLEXARNOLDICPMAPSALGORITHM_H
#define COMPLEXARNOLDICPMAPSALGORITHM_H 


#include "config.h"
#include "MPSObjects/CompletelyPositiveMap.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexUpperHessenbergMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class ComplexArnoldiCPMapsAlgorithm 
{

 protected:
  
  CompletelyPositiveMap * Map;
  double Accuracy;
  double GroundStateEnergy;
  
  // array where the vectors of the Krylov subspace
  ComplexMatrix * ArnoldiVectors;

  // maximum  number of iterations
  int MaximumNumberIteration;
  // current iteration index
  int Index;

  // garbage flag to avoid duplicating memory area
  GarbageFlag Flag;

  // Hessenberg matrix of the Hamiltonian in the Krylov subspace
  ComplexUpperHessenbergMatrix ReducedMatrix;
  // temporary matrix used to duplicated ReducedMatrix before diagonalize it
  ComplexUpperHessenbergMatrix TemporaryReducedMatrix;

  ComplexDiagonalMatrix ComplexDiagonalizedMatrix;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Arnoldi iteration
  Complex PreviousLastWantedEigenvalue;
  // value of the wanted eigenvalue at previous Arnoldi iteration
  Complex* ComplexPreviousWantedEigenvalues;
  // ground state 
  ComplexMatrix GroundState;

  // flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  bool StrongConvergenceFlag;
  //true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  bool HighEnergyFlag;
  // flag to compute the left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  bool LeftFlag;
  // sort the eigenvalues only with respect to their real part
  bool SortEigenvalueRealPartFlag;

  // array used to store temporary scalar products
  Complex* TemporaryCoefficients;

 public:

  // default constructor
  //
  // architecture = architecture to use for matrix operations
  // nbrEigenvalue = number of wanted eigenvalues
  // maxIter = an approximation of maximal number of iteration
  // highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
  // leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
  // strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
  ComplexArnoldiCPMapsAlgorithm (CompletelyPositiveMap * map, double accuracy, int nbrEigenvalue, int maxIter = 100, bool highEnergy = false, bool leftFlag = false, bool strongConvergence = false, bool sortRealFlag = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ComplexArnoldiCPMapsAlgorithm(const ComplexArnoldiCPMapsAlgorithm& algorithm);

  // destructor
  //
  ~ComplexArnoldiCPMapsAlgorithm();

  // initialize Arnoldi algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Arnoldi algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const ComplexMatrix & matrix);

  // get last produced vector
  //
  // return value = reference on last produced vector
  ComplexMatrix & GetGroundState();

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  ComplexMatrix * GetEigenstates(int nbrEigenstates);

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

 protected:

  // diagonalize tridiagonalized matrix and find ground state energy
  //
  virtual void Diagonalize();

  // restart the Arnoldi algorithm if needed
  //
  virtual void RestartAlgorithm();

};

#endif
