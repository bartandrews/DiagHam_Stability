////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//                   class of basic block Lanczos algorithm                   //
//                      (without re-orthogonalization )                       //
//                                                                            //
//                        last modification : 05/06/2008                      //
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


#include "LanczosAlgorithm/BasicBlockLanczosAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <iostream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
// blockSize = size of the block used for the block Lanczos algorithm
// maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
// lapackFlag = rely on LAPACK library to diagonalize the block matrix

BasicBlockLanczosAlgorithm::BasicBlockLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize, int maxIter,
						       bool diskFlag , bool resumeDiskFlag, bool strongConvergence, bool lapackFlag) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->BlockSize = blockSize;
  this->StrongConvergenceFlag = strongConvergence;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  if ((this->MaximumNumberIteration % this->BlockSize) != 0)
    this->MaximumNumberIteration = ((this->MaximumNumberIteration / this->BlockSize) + 1) * this->BlockSize;
  this->DiskFlag = diskFlag;
  this->ResumeDiskFlag = resumeDiskFlag;
  this->LanczosVectors = new RealVector [3 * this->BlockSize];
  this->TemporaryCoefficients = new double [3 * this->BlockSize];
  this->InitialStates = new RealVector [this->BlockSize];

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

BasicBlockLanczosAlgorithm::BasicBlockLanczosAlgorithm(const BasicBlockLanczosAlgorithm& algorithm) 
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
  this->DiskFlag = algorithm.DiskFlag;
  this->ResumeDiskFlag = algorithm.ResumeDiskFlag;
  this->LapackFlag = algorithm.LapackFlag;
}

// destructor
//

BasicBlockLanczosAlgorithm::~BasicBlockLanczosAlgorithm() 
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

void BasicBlockLanczosAlgorithm::InitializeLanczosAlgorithm() 
{
  if (this->ResumeDiskFlag == false)
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
      if (this->DiskFlag == false)
	{	  
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->InitialStates[j] = RealVector (this->LanczosVectors[j], true);
	}
      else
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      sprintf(TmpVectorName, "vector.%d", k);
	      this->LanczosVectors[k].WriteVector(TmpVectorName);
	    }
	  delete[] TmpVectorName;	  
	}
      this->Index = 0;
      this->ReducedMatrix.Resize(0, 0);
    }
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicBlockLanczosAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  if (this->ResumeDiskFlag == false)
    {
      int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
      this->LanczosVectors[0] = vector;
      double* TmpCoef = new double [this->NbrEigenvalue];
      for (int j = 1; j < this->NbrEigenvalue; ++j)
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
      if (this->DiskFlag == false)
	{	  
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->InitialStates[j] = RealVector (this->LanczosVectors[j], true);
	}
      else
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      sprintf(TmpVectorName, "vector.%d", k);
	      this->LanczosVectors[k].WriteVector(TmpVectorName);
	    }
	  delete[] TmpVectorName;	  
	}
      this->ReducedMatrix.Resize(0, 0);
    }
}

// initialize Lanczos algorithm with a set of given vectors
//
// vectors = array of vectors used as first step vectors
// nbrVectors = number of vectors in the array

void BasicBlockLanczosAlgorithm::InitializeLanczosAlgorithm(Vector* vectors, int nbrVectors)
{
  if (this->ResumeDiskFlag == false)
    {
      int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
      int MaxNbrVector = nbrVectors;
      if (MaxNbrVector > this->BlockSize)
	MaxNbrVector = this->BlockSize;
      for (int i = 0; i < MaxNbrVector; ++i)    
	this->LanczosVectors[i] = ((RealVector*) vectors)[i];
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
      if (this->DiskFlag == false)
	{	  
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->InitialStates[j] = RealVector (this->LanczosVectors[j], true);
	}
      else
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      sprintf(TmpVectorName, "vector.%d", k);
	      this->LanczosVectors[k].WriteVector(TmpVectorName);
	    }
	  delete[] TmpVectorName;	  
	}
      this->Index = 0;
      this->ReducedMatrix.Resize(0, 0);
    }
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& BasicBlockLanczosAlgorithm::GetGroundState()
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

Vector* BasicBlockLanczosAlgorithm::GetEigenstates(int nbrEigenstates)
{
  RealVector* Eigenstates = new RealVector [this->BlockSize];
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

  double* TmpCoefficents = new double [this->BlockSize];
  if (this->DiskFlag == false)
    {
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  Eigenstates[i].Copy(this->InitialStates[0], TmpEigenvector(0, i));
	  for (int j = 1; j < this->BlockSize; ++j)
	    TmpCoefficents[j - 1] = TmpEigenvector(j, i);	  
	  AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->InitialStates[1]), this->BlockSize - 1,  TmpCoefficents);
	  Operation.ApplyOperation(this->Architecture);
	  this->LanczosVectors[i].Copy(this->InitialStates[i]);
	}       
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
 	{
 	  for (int j = 0; j < this->BlockSize; ++j)
 	    TmpCoefficents[j] = TmpEigenvector(this->BlockSize + j, k);	  
 	  AddRealLinearCombinationOperation Operation5 (&(Eigenstates[k]), &(this->LanczosVectors[this->BlockSize]), this->BlockSize,  TmpCoefficents);
 	  Operation5.ApplyOperation(this->Architecture);
 	}       
      for (int i = 1; i < this->Index; ++i)
	{
	  this->ReducedMatrix.Resize((i + 1) * this->BlockSize, (i + 1) * this->BlockSize);
	  int NewVectorPosition = i * this->BlockSize;
	  MultipleVectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
								 &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      MultipleRealScalarProductOperation Operation (&(this->LanczosVectors[k + (2 * this->BlockSize)]), 
							    &(this->LanczosVectors[this->BlockSize]),   
							    k + 1, this->TemporaryCoefficients);
	      Operation.ApplyOperation(this->Architecture);
	      for (int j = 0; j <= k; ++j)
		{
		  this->ReducedMatrix(NewVectorPosition + k, NewVectorPosition + j) = this->TemporaryCoefficients[j];
		}
	    }
	  int Lim = (i - 1) * this->BlockSize;
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      for (int k = j; k < (2 * this->BlockSize); ++k)
		this->TemporaryCoefficients[k - j] = -this->ReducedMatrix(Lim + this->BlockSize + j, Lim + k);
	      AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[j + (2 * this->BlockSize)]), &(this->LanczosVectors[j]), 
							    2 * this->BlockSize - j,
							    this->TemporaryCoefficients);	  
	      Operation2.ApplyOperation(this->Architecture);
	    }
	  this->ReorthogonalizeVectors(&(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize, this->ReducedMatrix, 
				       NewVectorPosition - this->BlockSize, NewVectorPosition);  
 	  for (int k = 0; k < this->BlockSize; ++k)
 	    {
 	      for (int j = 0; j < this->BlockSize; ++j)
 		TmpCoefficents[j] = TmpEigenvector(((i + 1) * this->BlockSize) + j, k);	  
 	      AddRealLinearCombinationOperation Operation5 (&(Eigenstates[k]), &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize,  TmpCoefficents);
 	      Operation5.ApplyOperation(this->Architecture);
 	    }       	  
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      RealVector TmpVector = this->LanczosVectors[k];
	      this->LanczosVectors[k] = this->LanczosVectors[k + this->BlockSize];
	      this->LanczosVectors[k + this->BlockSize] = this->LanczosVectors[k + (2 * this->BlockSize)];
	      this->LanczosVectors[k + (2 * this->BlockSize)] = TmpVector;
	    }
	  cout << i << "/" << this->Index  << "           \r";
 	  cout.flush();
	}
    }
  else
    {
      char* TmpVectorName = new char [256];
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      sprintf(TmpVectorName, "vector.%d", j);
	      this->LanczosVectors[j].ReadVector(TmpVectorName);
	    }
	  Eigenstates[i].Copy(this->LanczosVectors[0], TmpEigenvector(0, i));
	  for (int j = 1; j < this->BlockSize; ++j)
	    TmpCoefficents[j - 1] = TmpEigenvector(j, i);	  
	  AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), this->BlockSize - 1,  TmpCoefficents);
	  Operation.ApplyOperation(this->Architecture);
	}       
      for (int i = 1; i < this->Index; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      sprintf(TmpVectorName, "vector.%d", ((i * this->BlockSize) + j));
	      this->LanczosVectors[j].ReadVector(TmpVectorName);
	    }
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      for (int j = 0; j < this->BlockSize; ++j)
		TmpCoefficents[j] = TmpEigenvector((i * this->BlockSize) + j, k);	  
	      AddRealLinearCombinationOperation Operation (&(Eigenstates[k]), this->LanczosVectors, this->BlockSize,  TmpCoefficents);
	      Operation.ApplyOperation(this->Architecture);
	    }
	  cout << i << "/" << this->Index << "           \r";
 	  cout.flush();
	}
      delete[] TmpVectorName;
    }
  for (int i = 0; i < this->BlockSize; ++i)
    Eigenstates[i] /= Eigenstates[i].Norm();
  cout << endl;
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicBlockLanczosAlgorithm::RunLanczosAlgorithm (int nbrIter) 
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
      if (this->DiskFlag == true)
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      sprintf(TmpVectorName, "vector.%d", (k + this->BlockSize));
	      this->LanczosVectors[k + this->BlockSize].WriteVector(TmpVectorName);
	    }
	  delete[] TmpVectorName;
	}
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

      if (this->ResumeDiskFlag == false)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      for (int k = j; k < (2 * this->BlockSize); ++k)
		this->TemporaryCoefficients[k - j] = -this->ReducedMatrix(Lim + this->BlockSize + j, Lim + k);
	      AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[j + (2 * this->BlockSize)]), &(this->LanczosVectors[j]), 
							    2 * this->BlockSize - j,
							    this->TemporaryCoefficients);	  
	      Operation2.ApplyOperation(this->Architecture);
	    }
	  this->ReorthogonalizeVectors(&(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize, this->ReducedMatrix, 
				       NewVectorPosition - this->BlockSize, NewVectorPosition);  
	}
      else
	{
	  this->ResumeDiskFlag = false;
	}
      if (this->DiskFlag == true)
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      sprintf(TmpVectorName, "vector.%d", (this->Index * this->BlockSize) + k);
	      this->LanczosVectors[(2 * this->BlockSize) + k].WriteVector(TmpVectorName);
	      RealVector TmpVector = this->LanczosVectors[k];
	      this->LanczosVectors[k] = this->LanczosVectors[k + this->BlockSize];
	      this->LanczosVectors[k + this->BlockSize] = this->LanczosVectors[k + (2 * this->BlockSize)];
	      this->LanczosVectors[k + (2 * this->BlockSize)] = TmpVector;
	      this->LanczosVectors[k] = RealVector();
	    }
	  delete[] TmpVectorName;
	  this->WriteState();
	}
      else
	{
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      RealVector TmpVector = this->LanczosVectors[k];
	      this->LanczosVectors[k] = this->LanczosVectors[k + this->BlockSize];
	      this->LanczosVectors[k + this->BlockSize] = this->LanczosVectors[k + (2 * this->BlockSize)];
	      this->LanczosVectors[k + (2 * this->BlockSize)] = TmpVector;
	    }
	}
      MultipleVectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
							     &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
      Operation2.ApplyOperation(this->Architecture);
      if (this->DiskFlag == true)
	{
	  char* TmpVectorName = new char [256];
	  for (int k = 0; k < this->BlockSize; ++k)
	    {
	      char* TmpVectorName = new char [256];
	      sprintf(TmpVectorName, "vector.%d", (k + ((this->Index - 1) * this->BlockSize)));
	      this->LanczosVectors[k].ReadVector(TmpVectorName);
	    }
	  delete[] TmpVectorName;
	}      
      for (int k = 0; k < this->BlockSize; ++k)
	{
	  MultipleRealScalarProductOperation Operation (&(this->LanczosVectors[k + (2 * this->BlockSize)]), 
							&(this->LanczosVectors[this->BlockSize]),   
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

void BasicBlockLanczosAlgorithm::Diagonalize () 
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

bool BasicBlockLanczosAlgorithm::TestConvergence ()
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

void BasicBlockLanczosAlgorithm::ReorthogonalizeVectors (RealVector* vectors, int nbrVectors, RealBandDiagonalSymmetricMatrix& matrix,
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

// read current Lanczos state from disk
//
// return value = true if no error occurs

bool BasicBlockLanczosAlgorithm::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->NbrEigenvalue);
  ReadLittleEndian(File, this->BlockSize);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->MaximumNumberIteration);  
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->ReducedMatrix.Resize(TmpDimension, TmpDimension);
  int TwiceBlockSize = 2 * this->BlockSize;
  int TmpMax = TmpDimension - TwiceBlockSize;
   for (int i = 0; i < TmpMax; ++i)    
    {    
      for (int j = 0; j < TwiceBlockSize; ++j)
	ReadLittleEndian(File, this->ReducedMatrix(i, i + j));
    }  
  for (int i = TmpMax; i < TmpDimension; ++i)    
    {    
      for (int j = TmpMax + 1; j < TmpDimension; ++j)
	ReadLittleEndian(File, this->ReducedMatrix(i, j));
    }  
  File.close();  
  char* TmpVectorName = new char [256];
  for (int k = 0; k < this->BlockSize; ++k)
    {
      sprintf(TmpVectorName, "vector.%d", ((this->Index * this->BlockSize) + k));
      this->LanczosVectors[k].ReadVector(TmpVectorName);
      sprintf(TmpVectorName, "vector.%d", (((this->Index + 1) * this->BlockSize) + k));
      this->LanczosVectors[k + this->BlockSize].ReadVector(TmpVectorName);
      sprintf(TmpVectorName, "vector.%d", (((this->Index + 2) * this->BlockSize) + k));
      this->LanczosVectors[k + (this->BlockSize * 2)].ReadVector(TmpVectorName);
    }
  delete[] TmpVectorName;
  return true;
}

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool BasicBlockLanczosAlgorithm::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->NbrEigenvalue);
  WriteLittleEndian(File, this->BlockSize);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->MaximumNumberIteration);  
  int TmpDimension = this->ReducedMatrix.GetNbrRow();
  WriteLittleEndian(File, TmpDimension);
  int TwiceBlockSize = 2 * this->BlockSize;
  int TmpMax = TmpDimension - TwiceBlockSize;
  for (int i = 0; i < TmpMax; ++i)    
    {    
      for (int j = 0; j < TwiceBlockSize; ++j)
	WriteLittleEndian(File, this->ReducedMatrix(i, i + j));
    }  
  for (int i = TmpMax; i < TmpDimension; ++i)    
    {    
      for (int j = TmpMax + 1; j < TmpDimension; ++j)
	WriteLittleEndian(File, this->ReducedMatrix(i, j));
    }  
  File.close();  
  return true;
}

