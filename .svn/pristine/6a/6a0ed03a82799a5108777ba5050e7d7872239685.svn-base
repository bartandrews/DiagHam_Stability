////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of basic Lanczos algorithm with complex vectors,         //
//                           ground state evaluation                          //
//                     and a projector over a set of vectors                  //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 11/06/2015                      //
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


#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateAndProjector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"

#include <stdlib.h>


// constructor
//
// nbrProjectors = dimension of the projector subspace
// projectorVectors = array that contains the vectors that spans the projector subspace
// projectorCoefficient = energy scale in front of the projector
// indexShiftFlag = true if the eigenstate indices have to be shifted
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration

ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::ComplexBasicLanczosAlgorithmWithGroundStateAndProjector(int nbrProjectors, ComplexVector* projectorVectors, double projectorCoefficient, bool indexShiftFlag,
														 AbstractArchitecture* architecture, int maxIter) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->V3 = ComplexVector();
  this->InitialState = ComplexVector();
  this->GroundStateFlag = false;
  this->NbrProjectors = nbrProjectors;
  this->ProjectorVectors = new ComplexVector [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    this->ProjectorVectors[i] = projectorVectors[i];
  this->ProjectorCoefficient = projectorCoefficient;
  this->IndexShiftFlag = indexShiftFlag;
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(maxIter, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(maxIter);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->Architecture = architecture;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::ComplexBasicLanczosAlgorithmWithGroundStateAndProjector(const ComplexBasicLanczosAlgorithmWithGroundStateAndProjector& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->V3 = algorithm.V3;
  this->InitialState = algorithm.InitialState;
  this->GroundStateFlag = algorithm.GroundStateFlag;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->NbrProjectors = algorithm.NbrProjectors;
  this->ProjectorVectors = new ComplexVector [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    this->ProjectorVectors[i] = algorithm.ProjectorVectors[i];
  this->ProjectorCoefficient = algorithm.ProjectorCoefficient;
  this->IndexShiftFlag = algorithm.IndexShiftFlag;
}

// destructor
//

ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::~ComplexBasicLanczosAlgorithmWithGroundStateAndProjector() 
{
}

// initialize Lanczos algorithm with a random vector
//

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = ComplexVector (Dimension);
  this->V2 = ComplexVector (Dimension);
  this->V3 = ComplexVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    {
      this->V1.Re(i) = (rand() - 32767) * 0.5;
      this->V1.Im(i) = (rand() - 32767) * 0.5;
    }
  this->V1 /= this->V1.Norm();
  Complex* TmpScalarProduct = new Complex[this->NbrProjectors];
  MultipleComplexScalarProductOperation Operation1 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation1.ApplyOperation(this->Architecture);	
  for (int i = 0; i < this->NbrProjectors; ++i)
    TmpScalarProduct[i] = -Conj(TmpScalarProduct[i]);
  AddComplexLinearCombinationOperation Operation2 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation2.ApplyOperation(this->Architecture);
  delete[] TmpScalarProduct;
  this->V1 /= this->V1.Norm();
  this->InitialState = ComplexVector (this->V1, true);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = vector;
  this->V2 = ComplexVector (Dimension);
  this->V3 = ComplexVector (Dimension);
  Complex* TmpScalarProduct = new Complex[this->NbrProjectors];
  MultipleComplexScalarProductOperation Operation1 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation1.ApplyOperation(this->Architecture);	
  for (int i = 0; i < this->NbrProjectors; ++i)
    TmpScalarProduct[i] = -Conj(TmpScalarProduct[i]);
  AddComplexLinearCombinationOperation Operation2 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation2.ApplyOperation(this->Architecture);
  delete[] TmpScalarProduct;
  this->V1 /= this->V1.Norm();
  this->InitialState = ComplexVector (this->V1, true);
  this->Index = 0;
  this->GroundStateFlag = false;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on last produced vector

Vector& ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::GetGroundState()
{
  if (this->GroundStateFlag == false)
    {
      RealVector TmpComponents (this->DiagonalizedMatrix.GetNbrRow());
      this->TridiagonalizedMatrix.Eigenvector(this->GroundStateEnergy, TmpComponents);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->InitialState, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
      this->AddProjectorContribution(this->InitialState, this->V3);
      this->V3.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(0), this->InitialState);
      this->V3 /= this->V3.Norm();
      this->V2 *= TmpComponents[this->DiagonalizedMatrix.GetNbrRow() - 1];
      this->V2.AddLinearCombination(TmpComponents[0], this->InitialState, TmpComponents[1], this->V3);
      this->V2.AddLinearCombination(TmpComponents[this->DiagonalizedMatrix.GetNbrRow() - 2], this->V1);
      ComplexVector TmpV (this->V2);
      this->V2 = this->InitialState; 
      this->InitialState = TmpV;
      int lim = this->DiagonalizedMatrix.GetNbrRow() - 3;
      cout << "Calculating vector: ";
      for (int i = 1; i < lim; ++i)
	{
	  VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &this->V3, &this->V1);
	  Operation2.ApplyOperation(this->Architecture);
	  this->AddProjectorContribution(this->InitialState, this->V3);
	  this->V1.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(i), this->V3, 
					-this->TridiagonalizedMatrix.UpperDiagonalElement(i - 1), this->V2);
	  this->V1 /= this->V1.Norm();
	  this->InitialState.AddLinearCombination(TmpComponents[i + 1], this->V1);
	  ComplexVector TmpV (this->V2);
	  this->V2 = this->V3;
	  this->V3 = this->V1;
	  this->V1 = TmpV;
	  cout << ".";
	  cout.flush();
	}
      cout << "done." << endl;
      this->InitialState /= this->InitialState.Norm();
//      cout << (this->InitialState * this->InitialState) << endl;
//      this->Hamiltonian->Multiply(this->InitialState, this->V3);
//      cout << "GroundStateEnergy =  " << (this->V3 * this->InitialState) << endl;
      this->GroundStateFlag = true;
    }
  return this->InitialState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V1, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->AddProjectorContribution(this->V1, this->V2);
       this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->V1 * this->V2).Re;
      this->V2.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index), 
				    this->V1);
      this->V2 /= this->V2.Norm(); 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &this->V2, &this->V3);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      this->V3.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1), V2, 
				    -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index), V1);
      this->V3 /= this->V3.Norm();
//      cout << this->V3 << endl;
      ComplexVector TmpV (this->V1);
      this->V1 = this->V2;
      this->V2 = this->V3;
      this->V3 = TmpV;
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V2, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->AddProjectorContribution(this->V2, this->V3);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
//      cout << this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) << endl;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
    }
  this->Diagonalize();
}



// add the projector contribution to the hamiltonian-vector multiplication
//
// initialVector = reference on the initial vector
// destinationVector = reference on the destination vector 

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjector::AddProjectorContribution(ComplexVector& initialVector, ComplexVector& destinationVector)
{
  Complex* TmpScalarProduct = new Complex[this->NbrProjectors];
  MultipleComplexScalarProductOperation Operation1 (&initialVector, this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation1.ApplyOperation(this->Architecture);	
  for (int i = 0; i < this->NbrProjectors; ++i)
    TmpScalarProduct[i] = this->ProjectorCoefficient * Conj(TmpScalarProduct[i]);
  AddComplexLinearCombinationOperation Operation2 (&destinationVector, this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
  Operation2.ApplyOperation(this->Architecture);
  delete[] TmpScalarProduct;
}

