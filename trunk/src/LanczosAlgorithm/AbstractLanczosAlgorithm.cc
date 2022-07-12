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

#include "Matrix/ComplexMatrix.h" //added by ba340
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"


// destructor
//

AbstractLanczosAlgorithm::~AbstractLanczosAlgorithm() 
{
}

// Set Hamiltonian on which Lanczos algorithm will be applied
//
// hamiltonian = Hamiltonian to use

void AbstractLanczosAlgorithm::SetHamiltonian (AbstractHamiltonian* hamiltonian)
{
  this->Hamiltonian = hamiltonian;
}

// get current diagonalized matrix
//
// return value = reference on current diagonalized matrix

RealTriDiagonalSymmetricMatrix& AbstractLanczosAlgorithm::GetDiagonalizedMatrix () 
{
  return this->DiagonalizedMatrix;
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues
void AbstractLanczosAlgorithm::GetEigenvalues (double*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new double [nbrEigenvalues];
  for (int i = 0; i < nbrEigenvalues; ++i)
    {
      eigenvalues[i] = this->DiagonalizedMatrix(i, i);
    }
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void AbstractLanczosAlgorithm::GetEigenvalues (Complex*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new Complex [nbrEigenvalues];
  for (int i = 0; i < nbrEigenvalues; ++i)
    {
      eigenvalues[i] = this->DiagonalizedMatrix(i, i);
    }
}

// get ground state energy
//
// return value = ground state energy

double AbstractLanczosAlgorithm::GetGroundStateEnergy()
{
  return this->GroundStateEnergy;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* AbstractLanczosAlgorithm::GetEigenstates(int nbrEigenstates)
{
  return 0;
}

// resume Lanczos algorithm from disk datas in current directory
//

void AbstractLanczosAlgorithm::ResumeLanczosAlgorithm()
{
  this->InitializeLanczosAlgorithm();
}
  
// initialize Lanczos algorithm with a set of given vectors
//
// vectors = array of vectors used as first step vectors
// nbrVectors = number of vectors in the array

void AbstractLanczosAlgorithm::InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors)
{
  this->InitializeLanczosAlgorithm(vectors[0]);
}

// force orthogonalization with respect to a set of vectors
//
// fileName = name of the file describing the set of vectors
// return value = true if no error occured

bool AbstractLanczosAlgorithm::ForceOrthogonalization(char* fileName)
{
  return false;
}

// diagonalize tridiagonalized matrix and find ground state energy
//

void AbstractLanczosAlgorithm::Diagonalize () 
{
  int Dimension = this->TridiagonalizedMatrix.GetNbrRow();
  this->DiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
  this->DiagonalizedMatrix.Diagonalize(50);
  this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(0);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (this->DiagonalizedMatrix.DiagonalElement(DiagPos) < this->GroundStateEnergy)
      this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(DiagPos);  
  return;
}

// test if convergence has been reached
//
// return value = true if convergence has been reached

bool AbstractLanczosAlgorithm::TestConvergence ()
{
  return true;
}

// optional shift of the eigenstate file name indices
//
// return value = index shift

int AbstractLanczosAlgorithm::EigenstateIndexShift()
{
  return 0;
}

// evaluate spectral response
//
// matrix = symmetric tridiagonal Hamiltonian matrix
// final_term = term at which continued fraction ends
// omega = input angular frequency
// term_start = term at which continued fraction starts (default=0)
// epsilon = small parameter to avoid poles (default=10^-10)
// return value = spectral response Green's function

Complex AbstractLanczosAlgorithm::EvaluateSpectralResponse(double omega, const double epsilon, int final_term, int term_start)
{     
  if (final_term<0) final_term = this->TridiagonalizedMatrix.GetNbrRow();
  Complex a, denom_fraction_numerator, current_denom; // for term == final_term
  
  current_denom = Complex(omega,epsilon)-this->TridiagonalizedMatrix(final_term-1,final_term-1);
  
  for (int n = final_term-2; n >term_start; --n) // n is index of Lanczos steps go backwards
  {
    a = this->TridiagonalizedMatrix(n, n);
    denom_fraction_numerator = this->TridiagonalizedMatrix(n-1,n)*this->TridiagonalizedMatrix(n-1,n);
    current_denom = Complex(omega,epsilon)- a - denom_fraction_numerator / current_denom;
  }
  
  return 1.0/current_denom;
}


// sample the spectral response and write to file
void AbstractLanczosAlgorithm::SampleSpectralResponse(std::ostream &Str, double omegaMin, double omegaMax, double epsilon, double omegaInterval, double spectralResolution)
{  
//double step=(omegaMax-omegaMin)/(nbrPoints);
//  double omega = omegaMin;
  Complex previous_response=0.0;
  double poleOmega, widthOmega=0;
  
  for (double omega = omegaMin; omega<omegaMax; omega+=omegaInterval)
  {
    // only print if the |difference| between adjacent points > spectralResolution (e.g. 1%)
    if (fabs((Norm(EvaluateSpectralResponse(omega, epsilon))-Norm(previous_response))/Norm(previous_response))>spectralResolution)
    {
      Str << omega <<" "<< Norm(EvaluateSpectralResponse(omega, epsilon)) << std::endl;
    }
    previous_response=EvaluateSpectralResponse(omega, epsilon);
  }
}

