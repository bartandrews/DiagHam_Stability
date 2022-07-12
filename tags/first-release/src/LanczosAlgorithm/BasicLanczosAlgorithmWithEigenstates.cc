////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of basic Lanczos algorithm                     //
//                          with access to eigenstates                        //
//                                                                            //
//                        last modification : 17/07/2001                      //
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


#include "LanczosAlgorithm/BasicLanczosAlgorithmWithEigenstates.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include <stdlib.h>


// default constructor
//
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration

BasicLanczosAlgorithmWithEigenstates::BasicLanczosAlgorithmWithEigenstates(AbstractArchitecture* architecture, int maxIter) 
{
  this->Index = 0;
  this->MaximumNumberIteration = maxIter;
  this->Hamiltonian = 0;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->Architecture = architecture;
  this->Flag.Initialize();
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

BasicLanczosAlgorithmWithEigenstates::BasicLanczosAlgorithmWithEigenstates(const BasicLanczosAlgorithmWithEigenstates& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
}

// destructor
//

BasicLanczosAlgorithmWithEigenstates::~BasicLanczosAlgorithmWithEigenstates() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
    }
}

// initialize Lanczos algorithm with a random vector
//

void BasicLanczosAlgorithmWithEigenstates::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = RealVector (Dimension);
  this->LanczosVectors[1] = RealVector (Dimension);
  this->LanczosVectors[2] = RealVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    this->LanczosVectors[0][i] = (rand() - 32767) * 0.5;
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicLanczosAlgorithmWithEigenstates::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = vector;
  this->LanczosVectors[1] = RealVector (Dimension);
  this->LanczosVectors[2] = RealVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get ground state
//
// return value = reference on ground state

Vector& BasicLanczosAlgorithmWithEigenstates::GetGroundState()
{
//  return this->LanczosVectors[Index + 1];
  RealVector TmpComponents (this->DiagonalizedMatrix.GetNbrRow());
  this->TridiagonalizedMatrix.Eigenvector(this->GroundStateEnergy, TmpComponents);
/*  RealVector TmpComponents2(TmpComponents, true);
  TmpComponents2 *= this->TridiagonalizedMatrix;
  TmpComponents2 /= this->GroundStateEnergy;*/
//  for (int i = 0; i < this->DiagonalizedMatrix.GetNbrRow(); i++)
//    cout << i << " : " << TmpComponents[i] << endl; 
  this->GroundState.Copy(this->LanczosVectors[0], TmpComponents[0]);
//  this->GroundState = RealVector(this->LanczosVectors[0].GetVectorDimension(), true);
  for (int i = 1; i < this->DiagonalizedMatrix.GetNbrRow(); i++)
    this->GroundState.AddLinearCombination (TmpComponents[i], this->LanczosVectors[i]);
//  cout << this->GroundState.Norm() << endl;
  this->GroundState /= this->GroundState.Norm();
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicLanczosAlgorithmWithEigenstates::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      this->Architecture->Multiply(this->Hamiltonian, this->LanczosVectors[0], this->LanczosVectors[1]);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->LanczosVectors[0] * 
							    this->LanczosVectors[1]);
      this->LanczosVectors[1].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index), 
						   this->LanczosVectors[0]);
      this->LanczosVectors[1] /= this->LanczosVectors[1].Norm(); 
      this->Architecture->Multiply(this->Hamiltonian, this->LanczosVectors[1], this->LanczosVectors[2]);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[0] * 
								       this->LanczosVectors[2]);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[1] * 
								      this->LanczosVectors[2]);
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      this->LanczosVectors[i].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index + 1), 
						   this->LanczosVectors[i - 1], 
						   -this->TridiagonalizedMatrix.
						   UpperDiagonalElement(this->Index), 
						   this->LanczosVectors[i - 2]);
      this->LanczosVectors[i] /= this->LanczosVectors[i].Norm();
      this->Index++;
      this->LanczosVectors[i + 1] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      this->Architecture->Multiply(this->Hamiltonian, this->LanczosVectors[i], this->LanczosVectors[i + 1]);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[i - 1] * 
								       this->LanczosVectors[i + 1]);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[i] * 
								      this->LanczosVectors[i + 1]);
    }
  this->Diagonalize();
}
  
