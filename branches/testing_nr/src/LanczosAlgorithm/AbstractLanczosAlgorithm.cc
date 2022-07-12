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


