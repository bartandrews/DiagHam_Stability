////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of basic Lanczos algorithm with complex vectors          //
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


#include "config.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithProjector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include <stdlib.h>


// default constructor
//
// nbrProjectors = dimension of the projector subspace
// projectorVectors = array that contains the vectors that spans the projector subspace
// projectorCoefficient = energy scale in front of the projector
// indexShiftFlag = true if the eigenstate indices have to be shifted
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration

ComplexBasicLanczosAlgorithmWithProjector::ComplexBasicLanczosAlgorithmWithProjector(int nbrProjectors, ComplexVector* projectorVectors, double projectorCoefficient, bool indexShiftFlag,
										     AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->NbrEigenvalue = nbrEigenvalue;
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
  this->PreviousLastWantedEigenvalue = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ComplexBasicLanczosAlgorithmWithProjector::ComplexBasicLanczosAlgorithmWithProjector(const ComplexBasicLanczosAlgorithmWithProjector& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->NbrProjectors = algorithm.NbrProjectors;
  this->ProjectorVectors = new ComplexVector [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    this->ProjectorVectors[i] = algorithm.ProjectorVectors[i];
  this->ProjectorCoefficient = algorithm.ProjectorCoefficient;
  this->IndexShiftFlag = algorithm.IndexShiftFlag;
}

// destructor
//

ComplexBasicLanczosAlgorithmWithProjector::~ComplexBasicLanczosAlgorithmWithProjector() 
{
}

// initialize Lanczos algorithm with a random vector
//

void ComplexBasicLanczosAlgorithmWithProjector::InitializeLanczosAlgorithm() 
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
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ComplexBasicLanczosAlgorithmWithProjector::InitializeLanczosAlgorithm(const Vector& vector) 
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
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on last produced vector

Vector& ComplexBasicLanczosAlgorithmWithProjector::GetGroundState()
{
  return this->V2;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ComplexBasicLanczosAlgorithmWithProjector::RunLanczosAlgorithm (int nbrIter) 
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
      this->AddProjectorContribution(this->V2, this->V3);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  ComplexVector* TmpVector = new ComplexVector[2];
  double* TmpCoefficient = new double[2];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      TmpVector[0] = this->V1;
      TmpVector[1] = this->V2;
      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
      AddComplexLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);
      this->V3 /= this->V3.Norm();
      ComplexVector TmpV (this->V1);
      this->V1 = this->V2;
      this->V2 = this->V3;
      this->V3 = TmpV;
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation3 (this->Hamiltonian, &this->V2, &this->V3);
      Operation3.ApplyOperation(this->Architecture);
      this->AddProjectorContribution(this->V2, this->V3);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
    }
  delete[] TmpVector;
  delete[] TmpCoefficient;
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
    }
}
  
// add the projector contribution to the hamiltonian-vector multiplication
//
// initialVector = reference on the initial vector
// destinationVector = reference on the destination vector 

void ComplexBasicLanczosAlgorithmWithProjector::AddProjectorContribution(ComplexVector& initialVector, ComplexVector& destinationVector)
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

