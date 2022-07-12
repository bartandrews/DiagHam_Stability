////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of basic block Arnoldi algorithm                   //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 05/12/2012                      //
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


#include "LanczosAlgorithm/BasicBlockArnoldiAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// blockSize = size of the block used for the block Lanczos algorithm
// maxIter = an approximation of maximal number of iteration
// highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
// leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

BasicBlockArnoldiAlgorithm::BasicBlockArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize, int maxIter, 
						       bool highEnergy, bool leftFlag, bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->ArnoldiVectors = new RealVector [this->MaximumNumberIteration];
  this->TemporaryCoefficients = new double [this->MaximumNumberIteration];
  this->BlockSize = blockSize;
//   if (this->BlockSize < 2)
//     this->BlockSize = 2;
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = RealMatrix(this->MaximumNumberIteration, this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix();
      this->ReducedMatrix = RealMatrix();
   }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->StrongConvergenceFlag = strongConvergence;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->ComplexPreviousWantedEigenvalues = new Complex [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->ComplexPreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
  this->HighEnergyFlag = highEnergy;
  this->LeftFlag = leftFlag;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

BasicBlockArnoldiAlgorithm::BasicBlockArnoldiAlgorithm(const BasicBlockArnoldiAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->ArnoldiVectors = new RealVector [this->MaximumNumberIteration];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->BlockSize = algorithm.BlockSize;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;
  this->ComplexDiagonalizedMatrix = algorithm.ComplexDiagonalizedMatrix;
  this->ReducedMatrix = algorithm.ReducedMatrix;
  this->ComplexPreviousWantedEigenvalues = new Complex [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->ComplexPreviousWantedEigenvalues[i] = 0.0;
  this->TemporaryCoefficients = algorithm.TemporaryCoefficients;
  this->HighEnergyFlag = algorithm.HighEnergyFlag;
  this->LeftFlag = algorithm.LeftFlag;
}

// destructor
//

BasicBlockArnoldiAlgorithm::~BasicBlockArnoldiAlgorithm() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ArnoldiVectors;
      delete[]  this->TemporaryCoefficients;
    }
  delete[] this->ComplexPreviousWantedEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void BasicBlockArnoldiAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  for (int i = 0; i < (3 * this->BlockSize); ++i)
    {
      this->ArnoldiVectors[i] = RealVector (Dimension);
    }
  double* TmpCoef = new double [this->BlockSize];
  for (int j = 0; j <  this->BlockSize; ++j)
    {
      for (int i = 0; i < Dimension; i++)
	{
	  this->ArnoldiVectors[j][i] = (rand() - 32767) * 0.5;
	}
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->ArnoldiVectors[j] * this->ArnoldiVectors[i];
      for (int i = 0; i < j; ++i)
	this->ArnoldiVectors[j].AddLinearCombination(-TmpCoef[i], this->ArnoldiVectors[i]);
      this->ArnoldiVectors[j] /= this->ArnoldiVectors[j].Norm();
    }
  delete[] TmpCoef;
  this->Index = 0;
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicBlockArnoldiAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  for (int i = 1; i < (3 * this->BlockSize); ++i)
    {
      this->ArnoldiVectors[i] = RealVector (Dimension);
    }
  this->ArnoldiVectors[0] = vector;
  double* TmpCoef = new double [this->BlockSize];
  for (int j = 1; j <  this->BlockSize; ++j)
    {
      for (int i = 0; i < Dimension; i++)
	{
	  this->ArnoldiVectors[j][i] = (rand() - 32767) * 0.5;
	}
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->ArnoldiVectors[j] * this->ArnoldiVectors[i];
      for (int i = 0; i < j; ++i)
	this->ArnoldiVectors[j].AddLinearCombination(-TmpCoef[i], this->ArnoldiVectors[i]);
      this->ArnoldiVectors[j] /= this->ArnoldiVectors[j].Norm();
    }
  this->Index = 0;
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& BasicBlockArnoldiAlgorithm::GetGroundState()
{
  this->GroundState = ComplexVector(this->ArnoldiVectors[0].GetLargeVectorDimension(), true);
  ComplexMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  TmpEigenvector.SetToIdentity();

  ComplexDiagonalMatrix  SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn(), true);
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector, this->LeftFlag);
  for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
    SortedDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "lapack is required for BasicBlockArnoldiAlgorithm" << endl;
#endif
  if (this->HighEnergyFlag == false)
    SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  else
    SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector);

  Complex* TmpCoefficents = new Complex [SortedDiagonalizedMatrix.GetNbrColumn()];
  for (int j = 0; j < SortedDiagonalizedMatrix.GetNbrColumn(); ++j)
    TmpCoefficents[j] = TmpEigenvector[0][j];
  AddComplexLinearCombinationOperation Operation (&(this->GroundState), this->ArnoldiVectors, 
						  SortedDiagonalizedMatrix.GetNbrColumn(),  TmpCoefficents);
  delete[] TmpCoefficents;
  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* BasicBlockArnoldiAlgorithm::GetEigenstates(int nbrEigenstates)
{
  ComplexVector* Eigenstates = new ComplexVector [nbrEigenstates];
  ComplexMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  TmpEigenvector.SetToIdentity();

  ComplexDiagonalMatrix  SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn(), true);
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector, this->LeftFlag);
  for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
    SortedDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "lapack is required for BasicBlockArnoldiAlgorithm" << endl;
#endif
  if (this->HighEnergyFlag == false)
    SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  else
    SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector);

  Complex* TmpCoefficents = new Complex [SortedDiagonalizedMatrix.GetNbrColumn()];
  for (int i = 0; i < nbrEigenstates; ++i)
    {
      for (int j = 0; j < SortedDiagonalizedMatrix.GetNbrColumn(); ++j)
	TmpCoefficents[j] = TmpEigenvector[i][j];
      Eigenstates[i] = ComplexVector(this->ArnoldiVectors[0].GetVectorDimension(), true);
      AddComplexLinearCombinationOperation Operation (&(Eigenstates[i]), this->ArnoldiVectors, 
						      SortedDiagonalizedMatrix.GetNbrColumn(),  TmpCoefficents);
      Operation.ApplyOperation(this->Architecture);
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicBlockArnoldiAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      if (nbrIter < 2)
	nbrIter = 2;
      Dimension = nbrIter * this->BlockSize;
      this->ReducedMatrix.Resize(Dimension, Dimension);

      MultipleVectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, this->ArnoldiVectors, &(this->ArnoldiVectors[this->BlockSize]), 
							     this->BlockSize);
      Operation1.ApplyOperation(this->Architecture);
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleRealScalarProductOperation Operation2 (&(this->ArnoldiVectors[i + this->BlockSize]), this->ArnoldiVectors,   
							 this->BlockSize, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      this->ReducedMatrix[i][j] = this->TemporaryCoefficients[j];
	    }
	}

      for (int i = 0; i < this->BlockSize; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->TemporaryCoefficients[j] = -this->ReducedMatrix(i, j);
	  AddRealLinearCombinationOperation Operation2 (&(this->ArnoldiVectors[this->BlockSize + i]), this->ArnoldiVectors, 
						       this->BlockSize, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	}

      this->ReorthogonalizeVectors(&(this->ArnoldiVectors[this->BlockSize]), this->BlockSize, this->ReducedMatrix, this->BlockSize, 0);

      MultipleVectorHamiltonianMultiplyOperation Operation3 (this->Hamiltonian, &(this->ArnoldiVectors[this->BlockSize]), 
							     &(this->ArnoldiVectors[2 * this->BlockSize]), this->BlockSize);
      Operation3.ApplyOperation(this->Architecture);
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleRealScalarProductOperation Operation2 (&(this->ArnoldiVectors[i + (2 * this->BlockSize)]), 
							 this->ArnoldiVectors,   
							 2 * this->BlockSize, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j < (2 * this->BlockSize); ++j)
	    {
	      this->ReducedMatrix[this->BlockSize + i][j] = this->TemporaryCoefficients[j];
	    }
	}
      nbrIter -= 2;
      this->Index = 2;
    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter * this->BlockSize;
      this->ReducedMatrix.Resize(Dimension, Dimension);
    }
  for (; nbrIter > 0; --nbrIter)
    {
      int NewVectorPosition = this->Index * this->BlockSize;
      for (int j = 0; j < this->BlockSize; ++j)
	{
	  for (int k = 0; k < NewVectorPosition; ++k)
	    {
	      this->ReducedMatrix.GetMatrixElement(k, NewVectorPosition - this->BlockSize + j, this->TemporaryCoefficients[k]);
	      this->TemporaryCoefficients[k] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation2 (&(this->ArnoldiVectors[NewVectorPosition + j]), this->ArnoldiVectors, 
							NewVectorPosition, this->TemporaryCoefficients);	  
	  Operation2.ApplyOperation(this->Architecture);
	}

      this->ReorthogonalizeVectors(&(this->ArnoldiVectors[NewVectorPosition]), this->BlockSize, this->ReducedMatrix, 
				   NewVectorPosition, NewVectorPosition - this->BlockSize);  

      for (int j = 0; j < this->BlockSize; ++j)
	this->ArnoldiVectors[NewVectorPosition + this->BlockSize + j] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->ArnoldiVectors[NewVectorPosition]), 
							    &(this->ArnoldiVectors[NewVectorPosition + this->BlockSize]), this->BlockSize);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < this->BlockSize; ++j)
	{
	  MultipleRealScalarProductOperation Operation3 (&(this->ArnoldiVectors[NewVectorPosition + this->BlockSize + j]), this->ArnoldiVectors,
							 NewVectorPosition + this->BlockSize, this->TemporaryCoefficients);
	  Operation3.ApplyOperation(this->Architecture);
	  for (int k = 0; k < (NewVectorPosition + this->BlockSize); ++k)
	    {
	      this->ReducedMatrix.SetMatrixElement(k, NewVectorPosition + j, this->TemporaryCoefficients[k]);
	    }
	}
      ++this->Index;
    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->ComplexPreviousWantedEigenvalues[i] = this->ComplexDiagonalizedMatrix[i];
      this->Diagonalize();
      if (this->HighEnergyFlag == false)
	this->ComplexDiagonalizedMatrix.SortMatrixUpOrder();
      else
	this->ComplexDiagonalizedMatrix.SortMatrixDownOrder();
    }
  else
    {
      this->Diagonalize();
      if (this->HighEnergyFlag == false)
	this->ComplexDiagonalizedMatrix.SortMatrixUpOrder();
      else
	this->ComplexDiagonalizedMatrix.SortMatrixDownOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->ComplexPreviousWantedEigenvalues[i] = 2.0 * this->ComplexDiagonalizedMatrix[i];
    }
}

  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool BasicBlockArnoldiAlgorithm::TestConvergence ()
{
  if (this->ReducedMatrix.GetNbrRow() > this->NbrEigenvalue)
    {
      cout << this->Index << " " << this->ComplexDiagonalizedMatrix[0] << " " << this->ComplexPreviousWantedEigenvalues[this->NbrEigenvalue - 1] << " " 
	   << Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1] - this->ComplexPreviousWantedEigenvalues[this->NbrEigenvalue - 1]) << endl;
      if (this->StrongConvergenceFlag == true)
	{
	  for (int i = this->NbrEigenvalue - 1; i >= 0; --i)
	    {
	      if (Norm(this->ComplexDiagonalizedMatrix[i] - this->ComplexPreviousWantedEigenvalues[i]) > 
		  (this->EigenvaluePrecision * Norm(this->ComplexDiagonalizedMatrix[i])))
		{ 
		  if (Norm(ComplexDiagonalizedMatrix[i]) > MACHINE_PRECISION)
		    return false;
		  else
		    if (Norm(this->ComplexPreviousWantedEigenvalues[i]) > MACHINE_PRECISION)
		      return false;
		}
	    }
	  return true;
	}
      else
	if (Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1] - this->ComplexPreviousWantedEigenvalues[this->NbrEigenvalue - 1]) < 
	    (this->EigenvaluePrecision * Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1])))
	  {
	    return true;
	  }
	else
	  return false;
    }
  return false;
}

// diagonalize tridiagonalized matrix and find ground state energy
//

void BasicBlockArnoldiAlgorithm::Diagonalize () 
{
  int Dimension = this->ReducedMatrix.GetNbrRow();
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
  if (this->Index < 5)
    cout << this->TemporaryReducedMatrix << endl;
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (this->TemporaryReducedMatrix.GetNbrColumn());
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag);
  this->ComplexDiagonalizedMatrix.Resize(this->TemporaryReducedMatrix.GetNbrColumn(), this->TemporaryReducedMatrix.GetNbrColumn());
  for (int i = 0; i < this->TemporaryReducedMatrix.GetNbrColumn(); ++i)
    this->ComplexDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "error, LAPACK is required for BasicBlockArnoldiAlgorithm" << endl;
#endif
  this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[0]);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (Norm(this->ComplexDiagonalizedMatrix[DiagPos]) < this->GroundStateEnergy)
      this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[DiagPos]);  
  return;
}

// reorthogonalize a set of vectors using Gram-Schmidt algorithm
//
// vectors = array of vectors to reorthogonalize
// nbrVectors = number of vectors to reorthogonalize
// matrix = matrix where transformation matrix has to be stored
// rowShift = shift to apply to matrix row index to reach the upper leftmost element
// columnShift = shift to apply to matrix column index to reach the upper leftmost element

void BasicBlockArnoldiAlgorithm::ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealMatrix& matrix,
							 int rowShift, int columnShift)
{
  double TmpNorm = vectors[0].Norm();
  matrix[columnShift][rowShift] = TmpNorm;
  vectors[0] /= TmpNorm;
  for (int i = 1; i < nbrVectors; ++i)
    {
      MultipleRealScalarProductOperation Operation (&(vectors[i]), 
						    vectors,   
						    i, this->TemporaryCoefficients);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < i; ++j)
	{
	  matrix(rowShift + j, columnShift + i) = this->TemporaryCoefficients[j];
	  this->TemporaryCoefficients[j] *= -1.0;
	}
      AddRealLinearCombinationOperation Operation2 (&(vectors[i]), vectors, 
						    i, this->TemporaryCoefficients);	  
      Operation2.ApplyOperation(this->Architecture);
      TmpNorm = vectors[i].Norm();      
      matrix[columnShift + i][rowShift + i] = TmpNorm;
      vectors[i] /= TmpNorm;      
    }
}

