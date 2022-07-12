////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of abstract Lanczos algorithm                   //
//                                                                            //
//                        last modification : 30/04/2001                      //
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


#ifndef ABSTRACTLANCZOSALGORITHM_H
#define ABSTRACTLANCZOSALGORITHM_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Matrix/ComplexMatrix.h" //added by ba340

class RealVector;
class ComplexVector;
class AbstractArchitecture;

class AbstractLanczosAlgorithm
{

 protected:

  RealTriDiagonalSymmetricMatrix TridiagonalizedMatrix;
  RealTriDiagonalSymmetricMatrix DiagonalizedMatrix;
  AbstractHamiltonian* Hamiltonian;

  double GroundStateEnergy;

  // architecture to use for matrix operations
  AbstractArchitecture* Architecture;

  // precision needed for eigenvalues
  double EigenvaluePrecision;

  // precision needed for eigenvectors
  double EigenvectorPrecision;

 public:

  // destructor
  //
  virtual ~AbstractLanczosAlgorithm();

  // Set Hamiltonian on which Lanczos algorithm will be applied
  //
  // hamiltonian = Hamiltonian to use
  virtual void SetHamiltonian (AbstractHamiltonian* hamiltonian);

  // initialize Lanczos algorithm with a random vector
  //
  virtual void InitializeLanczosAlgorithm() = 0;
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  virtual void InitializeLanczosAlgorithm(const Vector& vector) = 0;

  // initialize Lanczos algorithm with a set of given vectors
  //
  // vectors = array of vectors used as first step vectors
  // nbrVectors = number of vectors in the array
  virtual void InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors);

  // resume Lanczos algorithm from disk datas in current directory
  //
  virtual void ResumeLanczosAlgorithm();
  
  // force orthogonalization with respect to a set of vectors
  //
  // fileName = name of the file describing the set of vectors
  // return value = true if no error occured
  virtual bool ForceOrthogonalization(char* fileName);

  // get ground state
  //
  // return value = reference on ground state
  virtual Vector& GetGroundState() = 0;

  // get ground state energy
  //
  // return value = ground state energy
  virtual double GetGroundStateEnergy();

  // get the n first eigenstates
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  virtual Vector* GetEigenstates(int nbrEigenstates);

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

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  virtual void RunLanczosAlgorithm (int nbrIter) = 0;
  
  // get current diagonalized matrix
  //
  // return value = reference on current diagonalized matrix
  virtual RealTriDiagonalSymmetricMatrix& GetDiagonalizedMatrix ();

  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  virtual bool TestConvergence ();

  // set precision needed for eigenvalues
  //
  // precision = precision needed for eigenvalues
  virtual void SetEigenvaluePrecision (double precision);

  // set precision needed for eigenvectors
  //
  // precision = precision needed for eigenvectors
  virtual void SetEigenvectorPrecision (double precision);

  // optional shift of the eigenstate file name indices
  //
  // return value = index shift
  virtual int EigenstateIndexShift();
  
  // evaluate spectral response
  //
  // final_term = term at which continued fraction ends
  // omega = input angular frequency
  // term_start = term at which continued fraction starts (default=0)
  // epsilon = small parameter to avoid poles (default=10^-10)
  // return value = spectral response Green's function
  virtual Complex EvaluateSpectralResponse(double omega, const double epsilon=1E-10, int final_term=-1, int term_start=0);
  
  // sample the spectral response and write to file
  //
  // Str = output file
  // omegaMin = minimum omega
  // omegaMax = maximum omega
  // epsilon = epsilon small parameter to avoid poles
  // omegaInterval = omega step size
  // spectralResolution = only print points which differ from adjacent point by > spectralResoltion
  virtual void SampleSpectralResponse(std::ostream &Str, double omegaMin, double omegaMax, double epsilon, double omegaInterval, double spectralResolution);
  

 protected:
  
  // diagonalize tridiagonalized matrix and find ground state energy
  //
  virtual void Diagonalize ();

};

// set precision needed for eigenvalues
//
// precision = precision needed for eigenvalues

inline void AbstractLanczosAlgorithm::SetEigenvaluePrecision (double precision)
{
  this->EigenvaluePrecision = precision;
}

// set precision needed for eigenvectors
//
// precision = precision needed for eigenvectors

inline void AbstractLanczosAlgorithm::SetEigenvectorPrecision (double precision)
{
  this->EigenvectorPrecision = precision;
}

#endif
