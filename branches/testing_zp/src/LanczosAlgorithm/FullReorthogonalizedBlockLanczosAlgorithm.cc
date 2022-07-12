////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//           class of full reorthogonalized block Lanczos algorithm           //
//                (with full re-orthogonalization at each step)               //
//                                                                            //
//                        last modification : 10/03/2005                      //
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


#include "LanczosAlgorithm/FullReorthogonalizedBlockLanczosAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
// blockSize = size of the block used for the block Lanczos algorithm
// maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
// lapackFlag = rely on LAPACK library to diagonalize the block matrix

FullReorthogonalizedBlockLanczosAlgorithm::FullReorthogonalizedBlockLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize, int maxIter,
										     bool strongConvergence, bool lapackFlag) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->BlockSize = blockSize;
  this->StrongConvergenceFlag = strongConvergence;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  if ((this->MaximumNumberIteration % this->BlockSize) != 0)
    this->MaximumNumberIteration = ((this->MaximumNumberIteration / this->BlockSize) + 1) * this->BlockSize;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];
  this->TemporaryCoefficients = new double [this->MaximumNumberIteration];

  if (this->MaximumNumberIteration > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = RealBandDiagonalSymmetricMatrix(this->MaximumNumberIteration, this->BlockSize, true);
      this->TemporaryReducedMatrix = RealBandDiagonalSymmetricMatrix(this->MaximumNumberIteration, this->BlockSize, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ReducedMatrix = RealBandDiagonalSymmetricMatrix();
      this->TemporaryReducedMatrix = RealBandDiagonalSymmetricMatrix();
    }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->PreviousLastWantedEigenvalue = 0.0;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
  this->LapackFlag = lapackFlag;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

FullReorthogonalizedBlockLanczosAlgorithm::FullReorthogonalizedBlockLanczosAlgorithm(const FullReorthogonalizedBlockLanczosAlgorithm& algorithm) 
{
  this->LanczosVectors = algorithm.LanczosVectors;
  this->Index = algorithm.Index;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;

  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->BlockSize = algorithm.BlockSize;

  this->Hamiltonian = algorithm.Hamiltonian;
  this->Architecture = algorithm.Architecture;

  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = algorithm.PreviousWantedEigenvalues[i];

  this->TemporaryCoefficients = algorithm.TemporaryCoefficients;

  this->ReducedMatrix = algorithm.ReducedMatrix;
  this->TemporaryReducedMatrix = algorithm.TemporaryReducedMatrix;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->DiagonalizedMatrix = algorithm.DiagonalizedMatrix;

  this->Flag = algorithm.Flag;
  this->LapackFlag = algorithm.LapackFlag;
}

// destructor
//

FullReorthogonalizedBlockLanczosAlgorithm::~FullReorthogonalizedBlockLanczosAlgorithm() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
      delete[] this->TemporaryCoefficients;
    }
  delete[] this->PreviousWantedEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void FullReorthogonalizedBlockLanczosAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = RealVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    this->LanczosVectors[0][i] = (drand48() - 0.5) * 2.0;
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  double* TmpCoef = new double [this->NbrEigenvalue];
  for (int j = 1; j < this->BlockSize; ++j)
    {
      this->LanczosVectors[j] = RealVector (Dimension);
      for (int i = 0; i < Dimension; ++i)
	this->LanczosVectors[j][i] = (drand48() - 0.5) * 2.0;
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->LanczosVectors[j] * this->LanczosVectors[i];
      for (int i = 0; i < j; ++i)
	this->LanczosVectors[j].AddLinearCombination(-TmpCoef[i], this->LanczosVectors[i]);
      double TmpNorm = this->LanczosVectors[j].Norm();
      this->LanczosVectors[j] /= TmpNorm;
    }
  delete[] TmpCoef;
  this->Index = 0;
  this->ReducedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void FullReorthogonalizedBlockLanczosAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = vector;
  double* TmpCoef = new double [this->NbrEigenvalue];
  for (int j = 1; j < this->BlockSize; ++j)
    {
      this->LanczosVectors[j] = RealVector (Dimension);
      for (int i = 0; i < Dimension; ++i)
	this->LanczosVectors[j][i] = (drand48() - 0.5) * 2.0;
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->LanczosVectors[j] * this->LanczosVectors[i];
      for (int i = 0; i < j; ++i)
	  this->LanczosVectors[j].AddLinearCombination(-TmpCoef[i], this->LanczosVectors[i]);
      double TmpNorm = this->LanczosVectors[j].Norm();
      this->LanczosVectors[j] /= TmpNorm;
    }
  delete[] TmpCoef;
  this->Index = 0;
  this->ReducedMatrix.Resize(0, 0);
}

// initialize Lanczos algorithm with a set of given vectors
//
// vectors = array of vectors used as first step vectors
// nbrVectors = number of vectors in the array

void FullReorthogonalizedBlockLanczosAlgorithm::InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors)
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  int MaxNbrVector = nbrVectors;
  if (MaxNbrVector > this->BlockSize)
    MaxNbrVector = this->BlockSize;
  for (int i = 0; i < MaxNbrVector; ++i)    
    {
      this->LanczosVectors[i] = ((RealVector*) vectors)[i];
    }
  if (nbrVectors < this->BlockSize)
    {
      double* TmpCoef = new double [this->BlockSize];
      for (int j = nbrVectors; j < this->BlockSize; ++j)
	{
	  this->LanczosVectors[j] = RealVector (Dimension);
	  for (int i = 0; i < Dimension; ++i)
	    this->LanczosVectors[j][i] = (drand48() - 0.5) * 2.0;
	  for (int i = 0; i < j; ++i)
	    TmpCoef[i] = this->LanczosVectors[j] * this->LanczosVectors[i];
	  for (int i = 0; i < j; ++i)
	    this->LanczosVectors[j].AddLinearCombination(-TmpCoef[i], this->LanczosVectors[i]);
	  double TmpNorm = this->LanczosVectors[j].Norm();
	  this->LanczosVectors[j] /= TmpNorm;
	}
      delete[] TmpCoef;
    }
  this->Index = 0;
  this->ReducedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& FullReorthogonalizedBlockLanczosAlgorithm::GetGroundState()
{
  RealMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->ReducedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  if (this->LapackFlag == true)
    {
      RealDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
      this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
      for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
	SortedDiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
    }
  else
    {
#endif
      this->TemporaryReducedMatrix.Tridiagonalize(SortedDiagonalizedMatrix, 1e-7, TmpEigenvector);
      SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
#ifdef __LAPACK__
    }
#endif
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

  double* TmpCoefficents = new double [this->ReducedMatrix.GetNbrRow()];
  for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
    TmpCoefficents[j] = TmpEigenvector(j, 0);
  this->GroundState = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
  this->GroundState.Copy(this->LanczosVectors[0], TmpEigenvector(0, 0));
  AddRealLinearCombinationOperation Operation (&(this->GroundState), &(this->LanczosVectors[1]), this->ReducedMatrix.GetNbrRow() - 1, &(TmpCoefficents[1]));
  Operation.ApplyOperation(this->Architecture);
  this->GroundState /= this->GroundState.Norm();
  delete[] TmpCoefficents;

  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* FullReorthogonalizedBlockLanczosAlgorithm::GetEigenstates(int nbrEigenstates)
{
  RealVector* Eigenstates = new RealVector [nbrEigenstates];
  RealMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->ReducedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  if (this->LapackFlag == true)
    {
      RealDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
      this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
      for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
	SortedDiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
    }
  else
    {
#endif
      this->TemporaryReducedMatrix.Tridiagonalize(SortedDiagonalizedMatrix, 1e-7, TmpEigenvector);
      SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
#ifdef __LAPACK__
    }
#endif
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  double* TmpCoefficents = new double [this->ReducedMatrix.GetNbrRow()];
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
	TmpCoefficents[j] = TmpEigenvector(j, i);
      Eigenstates[i] = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
      Eigenstates[i].Copy(this->LanczosVectors[0], TmpEigenvector(0, i));
      AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), this->ReducedMatrix.GetNbrRow() - 1, &(TmpCoefficents[1]));
      Operation.ApplyOperation(this->Architecture);
      Eigenstates[i] /= Eigenstates[i].Norm();
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void FullReorthogonalizedBlockLanczosAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      if (nbrIter < 2)
	nbrIter = 2;
      Dimension = nbrIter * this->BlockSize;
      this->ReducedMatrix.Resize(Dimension, Dimension);

      for (int k = 0; k < this->BlockSize; ++k)
	this->LanczosVectors[k + this->BlockSize] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, this->LanczosVectors, &(this->LanczosVectors[this->BlockSize]), this->BlockSize);
      Operation.ApplyOperation(this->Architecture);
 
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleRealScalarProductOperation Operation2 (&(this->LanczosVectors[i + this->BlockSize]), this->LanczosVectors,   
							i + 1, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= i; ++j)
	    {
	      this->ReducedMatrix(i, j) = this->TemporaryCoefficients[j];
	    }
	}
      
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->TemporaryCoefficients[j] = -this->ReducedMatrix(i, j);
	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[this->BlockSize + i]), this->LanczosVectors, 
						       this->BlockSize, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	}

      this->ReorthogonalizeVectors(&(this->LanczosVectors[this->BlockSize]), this->BlockSize, this->ReducedMatrix, 0, this->BlockSize);
      
      for (int k = 0; k < this->BlockSize; ++k)
	this->LanczosVectors[k + (2 * this->BlockSize)] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation3 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
							     &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
      Operation3.ApplyOperation(this->Architecture);
 
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleRealScalarProductOperation Operation2 (&(this->LanczosVectors[i + (2 * this->BlockSize)]), 
							&(this->LanczosVectors[this->BlockSize]),   
							i + 1, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= i; ++j)
	    {
	      this->ReducedMatrix(this->BlockSize + i, this->BlockSize + j) = this->TemporaryCoefficients[j];
	    }
	}
      nbrIter -= 2;
      this->Index = 2;
    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + (nbrIter * this->BlockSize);
      this->ReducedMatrix.Resize(Dimension, Dimension);
    }
  for (; nbrIter > 0; --nbrIter)
    {
      int NewVectorPosition = this->Index * this->BlockSize;
      int Lim = (this->Index - 2) * this->BlockSize;

      for (int j = 0; j < this->BlockSize; ++j)
	{
	  for (int k = j; k < (2 * this->BlockSize); ++k)
	    this->TemporaryCoefficients[k - j] = -this->ReducedMatrix(Lim + this->BlockSize + j, Lim + k);
	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[j + NewVectorPosition]), &(this->LanczosVectors[Lim + j]), 
							2 * this->BlockSize - j,
							this->TemporaryCoefficients);	  
	  Operation2.ApplyOperation(this->Architecture);
	}

      if (this->Index > 2)
	{
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[k + NewVectorPosition]), this->LanczosVectors, Lim, 
							     this->TemporaryCoefficients);
	      Operation4.ApplyOperation(this->Architecture);
	      for (int j = 0; j < Lim; j++)
		{
		  this->TemporaryCoefficients[j] *= -1.0;
		}
	      AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[k + NewVectorPosition]), this->LanczosVectors, Lim,
							    this->TemporaryCoefficients);
	      Operation2.ApplyOperation(this->Architecture);
	    }
	}

      this->ReorthogonalizeVectors(&(this->LanczosVectors[NewVectorPosition]), this->BlockSize, this->ReducedMatrix, 
				   NewVectorPosition - this->BlockSize, NewVectorPosition);  


      for (int k = 0; k < this->BlockSize; ++k)
	this->LanczosVectors[k + this->BlockSize + NewVectorPosition] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[NewVectorPosition]), 
							     &(this->LanczosVectors[this->BlockSize + NewVectorPosition]), this->BlockSize);
      Operation2.ApplyOperation(this->Architecture);
 
      for (int k = 0; k < this->BlockSize; ++k)
	{
	  MultipleRealScalarProductOperation Operation (&(this->LanczosVectors[k + NewVectorPosition + this->BlockSize]), 
							&(this->LanczosVectors[NewVectorPosition]),   
							k + 1, this->TemporaryCoefficients);
	  Operation.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= k; ++j)
	    {
	      this->ReducedMatrix(NewVectorPosition + k, NewVectorPosition + j) = this->TemporaryCoefficients[j];
	    }
	}
      ++this->Index;
    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
       for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->PreviousWantedEigenvalues[i] = this->DiagonalizedMatrix.DiagonalElement(i);
     this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->PreviousWantedEigenvalues[i] = 2.0 * this->DiagonalizedMatrix.DiagonalElement(i);
    }
}

  
// diagonalize tridiagonalized matrix and find ground state energy
//

void FullReorthogonalizedBlockLanczosAlgorithm::Diagonalize () 
{
  int Dimension = this->ReducedMatrix.GetNbrRow();
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  if (this->LapackFlag == true)
    {
      RealDiagonalMatrix TmpDiag (this->TemporaryReducedMatrix.GetNbrColumn());
      this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag);
      this->DiagonalizedMatrix.Resize(this->TemporaryReducedMatrix.GetNbrColumn(), this->TemporaryReducedMatrix.GetNbrColumn());
      for (int i = 0; i < this->TemporaryReducedMatrix.GetNbrColumn(); ++i)
	this->DiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
    }
  else
    {
#endif
      this->TemporaryReducedMatrix.Tridiagonalize(this->DiagonalizedMatrix, 1e-7);
      this->DiagonalizedMatrix.Diagonalize();
#ifdef __LAPACK__
    }
#endif
  this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(0);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (this->DiagonalizedMatrix.DiagonalElement(DiagPos) < this->GroundStateEnergy)
      this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(DiagPos);  
  return;
}

// test if convergence has been reached
//
// return value = true if convergence has been reached

bool FullReorthogonalizedBlockLanczosAlgorithm::TestConvergence ()
{
  if (this->DiagonalizedMatrix.GetNbrRow() > this->NbrEigenvalue)
    {
      if (this->StrongConvergenceFlag == true)
	{
	  for (int i = this->NbrEigenvalue - 1; i >= 0; --i)
	    {
	      if (fabs(this->DiagonalizedMatrix.DiagonalElement(i) - this->PreviousWantedEigenvalues[i]) > 
		  (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(i))))
		{
		  return false;
		}
	    }
	  return true;
	}
      else
	{
	  if (fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
	      (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1))))
	    {
	      return true;
	    }
	  else
	    {
	      return false;
	    }
	}
    }
  return false;
}

// reorthogonalize a set of vectors using Gram-Schmidt algorithm
//
// vectors = array of vectors to reorthogonalize
// nbrVectors = number of vectors to reorthogonalize
// matrix = matrix where transformation matrix has to be stored
// rowShift = shift to apply to matrix row index to reach the upper leftmost element
// columnShift = shift to apply to matrix column index to reach the upper leftmost element

void FullReorthogonalizedBlockLanczosAlgorithm::ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealBandDiagonalSymmetricMatrix& matrix,
                                                                        int rowShift, int columnShift)
{
  double TmpNorm = vectors[0].Norm();
  matrix(rowShift, columnShift) = TmpNorm;
  vectors[0] /= TmpNorm;
  for (int i = 1; i < nbrVectors; ++i)
    {
      MultipleRealScalarProductOperation Operation (&(vectors[i]), 
						    vectors,   
						    i, this->TemporaryCoefficients);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < i; ++j)
	{
	  matrix(rowShift + i, columnShift + j) = this->TemporaryCoefficients[j];
	  this->TemporaryCoefficients[j] *= -1.0;
	}
      AddRealLinearCombinationOperation Operation2 (&(vectors[i]), vectors, 
						    i, this->TemporaryCoefficients);	  
      Operation2.ApplyOperation(this->Architecture);
      TmpNorm = vectors[i].Norm();      
      matrix(rowShift + i, columnShift + i) = TmpNorm;
      vectors[i] /= TmpNorm;      
    }
}
