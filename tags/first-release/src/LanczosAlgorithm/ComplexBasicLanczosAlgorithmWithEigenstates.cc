////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of basic Lanczos algorithm with complex vectors          //
//                           and access to eigenstates                        //
//                                                                            //
//                        last modification : 26/03/2002                      //
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


#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"

#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration

ComplexBasicLanczosAlgorithmWithEigenstates::ComplexBasicLanczosAlgorithmWithEigenstates(AbstractArchitecture* architecture, int maxIter) 
{
  this->Index = 0;
  this->MaximumNumberIteration = maxIter;
  this->Hamiltonian = 0;
  this->LanczosVectors = new ComplexVector [this->MaximumNumberIteration];
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
  this->Flag.Initialize();
  this->Architecture = architecture;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ComplexBasicLanczosAlgorithmWithEigenstates::ComplexBasicLanczosAlgorithmWithEigenstates(const ComplexBasicLanczosAlgorithmWithEigenstates& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->LanczosVectors = new ComplexVector [this->MaximumNumberIteration];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
}

// destructor
//

ComplexBasicLanczosAlgorithmWithEigenstates::~ComplexBasicLanczosAlgorithmWithEigenstates() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
    }
}

// initialize Lanczos algorithm with a random vector
//

void ComplexBasicLanczosAlgorithmWithEigenstates::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = ComplexVector (Dimension);
  this->LanczosVectors[1] = ComplexVector (Dimension);
  this->LanczosVectors[2] = ComplexVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    {
      this->LanczosVectors[0].Re(i) = (rand() - 32767) * 0.5;
      this->LanczosVectors[0].Im(i) = (rand() - 32767) * 0.5;
    }
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ComplexBasicLanczosAlgorithmWithEigenstates::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = vector;
  this->LanczosVectors[1] = ComplexVector (Dimension);
  this->LanczosVectors[2] = ComplexVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get ground state
//
// return value = reference on ground state

Vector& ComplexBasicLanczosAlgorithmWithEigenstates::GetGroundState()
{
  RealVector TmpComponents (this->DiagonalizedMatrix.GetNbrRow());
  this->TridiagonalizedMatrix.Eigenvector(this->GroundStateEnergy, TmpComponents);
  this->GroundState.Copy(this->LanczosVectors[0], TmpComponents[0]);
  for (int i = 1; i < this->DiagonalizedMatrix.GetNbrRow(); i++)
    {
      this->GroundState.AddLinearCombination (TmpComponents[i], this->LanczosVectors[i]);
    }
  this->GroundState /= this->GroundState.Norm();
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ComplexBasicLanczosAlgorithmWithEigenstates::RunLanczosAlgorithm (int nbrIter) 
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
							    this->LanczosVectors[1]).Re;
      this->LanczosVectors[1].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index), 
						   this->LanczosVectors[0]);
      this->LanczosVectors[1] /= this->LanczosVectors[1].Norm(); 
//      this->Hamiltonian->Multiply(this->LanczosVectors[1], this->LanczosVectors[2]);
      this->Architecture->Multiply(this->Hamiltonian, this->LanczosVectors[1], this->LanczosVectors[2]);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[0] * 
								       this->LanczosVectors[2]).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[1] * 
								      this->LanczosVectors[2]).Re;
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  if (Dimension > this->MaximumNumberIteration)
    {
      cout << "warning: too much iterations" << endl;
      return;
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
      this->LanczosVectors[i + 1] = ComplexVector(this->Hamiltonian->GetHilbertSpaceDimension());
      this->Architecture->Multiply(this->Hamiltonian, this->LanczosVectors[i], this->LanczosVectors[i + 1]);
//      this->Hamiltonian->Multiply(this->LanczosVectors[i], this->LanczosVectors[i + 1]);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[i - 1] * 
								       this->LanczosVectors[i + 1]).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[i] * 
								      this->LanczosVectors[i + 1]).Re;
    }
  this->Diagonalize();
}
  
