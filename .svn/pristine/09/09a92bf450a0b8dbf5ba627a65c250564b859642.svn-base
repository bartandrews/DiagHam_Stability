////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//         class of full reorthogonalized complex block Lanczos algorithm     //
//                (with full re-orthogonalization at each step)               //
//                                                                            //
//                        last modification : 05/01/2012                      //
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


#include "LanczosAlgorithm/FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;
using std::ios;


// basic constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
// blockSize = size of the block used for the block Lanczos algorithm
// maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 
// lapackFlag = rely on LAPACK library to diagonalize the block matrix

FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize, int maxIter,
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
  this->LanczosVectors = new ComplexVector [3 * this->BlockSize];
  this->TemporaryCoefficients = new Complex [this->MaximumNumberIteration];

  if (this->MaximumNumberIteration > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = HermitianMatrix(this->MaximumNumberIteration, true);
      this->TemporaryReducedMatrix = HermitianMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ReducedMatrix = HermitianMatrix();
      this->TemporaryReducedMatrix = HermitianMatrix();
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
  this->ResumeDiskFlag = false;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage(const FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage& algorithm) 
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
  this->ResumeDiskFlag = algorithm.ResumeDiskFlag;
}

// destructor
//

FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::~FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage() 
{
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::GetGroundState()
{
  ComplexMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
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
      cout <<  "error, Lapack is required" << endl;
#ifdef __LAPACK__
    }
#endif
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

  Complex* TmpCoefficents = new Complex [this->ReducedMatrix.GetNbrRow()];
  for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
    TmpCoefficents[j] = TmpEigenvector(j, 0);
  this->GroundState = ComplexVector (this->Hamiltonian->GetHilbertSpaceDimension());
  this->ReadLanczosVectors(this->LanczosVectors, 0, 0, this->BlockSize);
  this->GroundState.Copy(this->LanczosVectors[0], TmpCoefficents[0]);
  AddComplexLinearCombinationOperation Operation (&(this->GroundState), &(this->LanczosVectors[1]), this->BlockSize - 1, &(TmpCoefficents[1]));
  Operation.ApplyOperation(this->Architecture);
  for (int j = this->BlockSize; j < this->ReducedMatrix.GetNbrRow(); j += this->BlockSize)
    {
      this->ReadLanczosVectors(this->LanczosVectors, j, 0, this->BlockSize);
      AddComplexLinearCombinationOperation Operation2 (&(this->GroundState), this->LanczosVectors, this->BlockSize, TmpCoefficents);
      Operation2.ApplyOperation(this->Architecture);
    }
  this->GroundState /= this->GroundState.Norm();
  delete[] TmpCoefficents;

  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::GetEigenstates(int nbrEigenstates)
{
  ComplexVector* Eigenstates = new ComplexVector [nbrEigenstates];
  ComplexMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  TmpEigenvector.SetToIdentity();

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  if (this->LapackFlag == true)
    {
      RealDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
      this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
      for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
	{
	  SortedDiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
	}
    }
  else
    {
#endif
      RealDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
      this->TemporaryReducedMatrix.Diagonalize(TmpDiag, TmpEigenvector);
      for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
	{
	  SortedDiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
	}
#ifdef __LAPACK__
    }
#endif
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

  Complex** TmpCoefficents = new Complex* [this->ReducedMatrix.GetNbrRow()];
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      TmpCoefficents[i] = new Complex [this->ReducedMatrix.GetNbrRow()];
      for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
	TmpEigenvector.GetMatrixElement(j, i, TmpCoefficents[i][j]);
      Eigenstates[i] = ComplexVector (this->Hamiltonian->GetHilbertSpaceDimension());
    }
  this->ReadLanczosVectors(this->LanczosVectors, 0, 0, this->BlockSize);
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      Eigenstates[i].Copy(this->LanczosVectors[0], TmpCoefficents[i][0]);
      AddComplexLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), this->BlockSize - 1, &(TmpCoefficents[i][1]));
      Operation.ApplyOperation(this->Architecture);
    }
  for (int j = this->BlockSize; j < this->ReducedMatrix.GetNbrRow(); j += this->BlockSize)
    {
      this->ReadLanczosVectors(this->LanczosVectors, j, 0, this->BlockSize);
      for (int i = 0; i < nbrEigenstates; ++i) 
	{
	  AddComplexLinearCombinationOperation Operation2 (&(Eigenstates[i]), this->LanczosVectors, this->BlockSize, &(TmpCoefficents[i][j]));
	  Operation2.ApplyOperation(this->Architecture);
	}
    }
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      Eigenstates[i] /= Eigenstates[i].Norm();
      delete[] TmpCoefficents[i];
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// resume Lanczos algorithm from disk datas in current directory
//

void FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::ResumeLanczosAlgorithm()
{
  this->ReadState();
  this->ReadLanczosVectors(this->LanczosVectors, (this->Index - 2) * this->BlockSize, 0, 3 * this->BlockSize);
  if (this->Hamiltonian != NULL)
    {
      for (int i = 0; i < (3 * this->BlockSize); ++i)
	{
	  if (this->Hamiltonian->GetHilbertSpaceDimension() != this->LanczosVectors[i].GetVectorDimension())
	    {
	      cout << "Hamiltonian does not match dimension of stored vectors in resuming"<<endl;
	      exit(-1);
	    }
	}
    }
  this->ResumeDiskFlag = true;
}
  
// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  Complex Tmp;
  if (this->Index == 0)
    {
      if (nbrIter < 2)
	nbrIter = 2;
      Dimension = nbrIter * this->BlockSize;

      this->ReducedMatrix.Resize(Dimension, Dimension);

      for (int k = 0; k < this->BlockSize; ++k)
	this->LanczosVectors[k + this->BlockSize] = ComplexVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, this->LanczosVectors, &(this->LanczosVectors[this->BlockSize]), this->BlockSize);
      Operation.ApplyOperation(this->Architecture);
      
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleComplexScalarProductOperation Operation2 (&(this->LanczosVectors[i + this->BlockSize]), this->LanczosVectors,   
							    i + 1, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= i; ++j)
	    {
	      this->ReducedMatrix.SetMatrixElement(j, i, this->TemporaryCoefficients[j]);
	    }
	}


      for (int i = 0; i < this->BlockSize; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      this->ReducedMatrix.GetMatrixElement(i, j, this->TemporaryCoefficients[j]);
	      this->TemporaryCoefficients[j] *= -1.0;
	    }
	  AddComplexLinearCombinationOperation Operation2 (&(this->LanczosVectors[this->BlockSize + i]), this->LanczosVectors, 
							   this->BlockSize, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	}

      this->ReorthogonalizeVectors(&(this->LanczosVectors[this->BlockSize]), this->BlockSize, this->ReducedMatrix, 0, this->BlockSize);
      this->WriteLanczosVectors(this->LanczosVectors, this->BlockSize, this->BlockSize, this->BlockSize);
      
      for (int k = 0; k < this->BlockSize; ++k)
	this->LanczosVectors[k + (2 * this->BlockSize)] = ComplexVector(this->Hamiltonian->GetHilbertSpaceDimension());
      MultipleVectorHamiltonianMultiplyOperation Operation3 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
							     &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
      Operation3.ApplyOperation(this->Architecture);
 
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleComplexScalarProductOperation Operation2 (&(this->LanczosVectors[i + (2 * this->BlockSize)]), 
							    &(this->LanczosVectors[this->BlockSize]),   
							    i + 1, this->TemporaryCoefficients);
	  Operation2.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= i; ++j)
	    {
	      this->ReducedMatrix.SetMatrixElement(this->BlockSize + j, this->BlockSize + i, this->TemporaryCoefficients[j]);
	    }
	}

      nbrIter -= 2;
      this->Index = 2;
   }
  else
    {
      if (this->ResumeDiskFlag == false)
	{
	  Dimension = this->ReducedMatrix.GetNbrRow() + (nbrIter * this->BlockSize);
	  this->ReducedMatrix.Resize(Dimension, Dimension);
	}
    }
  for (; nbrIter > 0; --nbrIter)
    {
      int NewVectorPosition = this->Index * this->BlockSize;
      int VectorPositionShift = 2 * this->BlockSize;
      int Lim = (this->Index - 2) * this->BlockSize;

      if (this->ResumeDiskFlag == false)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    {
	      for (int k = j; k < (2 * this->BlockSize); ++k)
		{
		  this->ReducedMatrix.GetMatrixElement(Lim + this->BlockSize + j, Lim + k, this->TemporaryCoefficients[k - j]);
		  this->TemporaryCoefficients[k - j] *= -1.0;
		}
	      AddComplexLinearCombinationOperation Operation2 (&(this->LanczosVectors[j + VectorPositionShift]), &(this->LanczosVectors[j]), 
							       2 * this->BlockSize - j, this->TemporaryCoefficients);	  
	      Operation2.ApplyOperation(this->Architecture);
	    }
	  
	  for (int l = 0; l < Lim; l += this->BlockSize)
	    {
	      this->ReadLanczosVectors(this->LanczosVectors, l, 0, this->BlockSize);
	      for (int k = 0; k < this->BlockSize; ++k)
		{
		  MultipleComplexScalarProductOperation Operation4 (&(this->LanczosVectors[k + VectorPositionShift]), this->LanczosVectors, this->BlockSize, 
								    this->TemporaryCoefficients);
		  Operation4.ApplyOperation(this->Architecture);
		  for (int j = 0; j < Lim; j++)
		    {
		      this->TemporaryCoefficients[j] = Conj(this->TemporaryCoefficients[j]);
		      this->TemporaryCoefficients[j] *= -1.0;
		    }
		  AddComplexLinearCombinationOperation Operation2 (&(this->LanczosVectors[k + VectorPositionShift]), this->LanczosVectors, this->BlockSize,
								   this->TemporaryCoefficients);
		  Operation2.ApplyOperation(this->Architecture);
		}
	    }
	  
	  this->ReorthogonalizeVectors(&(this->LanczosVectors[VectorPositionShift]), this->BlockSize, this->ReducedMatrix, 
				       NewVectorPosition - this->BlockSize, NewVectorPosition);  
	  
	  this->WriteLanczosVectors(this->LanczosVectors, NewVectorPosition, 2 * this->BlockSize, this->BlockSize);
	  this->WriteState();
	}
      else
	{
	  this->ResumeDiskFlag = false;
	}
      for (int k = 0; k < this->BlockSize; ++k)
	{
	  ComplexVector TmpVector = this->LanczosVectors[k];
	  this->LanczosVectors[k] = this->LanczosVectors[k + this->BlockSize];
	  this->LanczosVectors[k + this->BlockSize] = this->LanczosVectors[k + (2 * this->BlockSize)];
	  this->LanczosVectors[k + (2 * this->BlockSize)] = TmpVector;
	}

      MultipleVectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
							     &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
      Operation2.ApplyOperation(this->Architecture);
 
      for (int k = 0; k < this->BlockSize; ++k)
	{
	  MultipleComplexScalarProductOperation Operation (&(this->LanczosVectors[k + (2* this->BlockSize)]), 
							   &(this->LanczosVectors[this->BlockSize]),   
							   k + 1, this->TemporaryCoefficients);
	  Operation.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= k; ++j)
	    {
	      this->TemporaryCoefficients[j] = Conj(this->TemporaryCoefficients[j]);
	      this->ReducedMatrix.SetMatrixElement(NewVectorPosition + k, NewVectorPosition + j, this->TemporaryCoefficients[j]);
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

  
// reorthogonalize a set of vectors using Gram-Schmidt algorithm
//
// vectors = array of vectors to reorthogonalize
// nbrVectors = number of vectors to reorthogonalize
// matrix = matrix where transformation matrix has to be stored
// rowShift = shift to apply to matrix row index to reach the upper leftmost element
// columnShift = shift to apply to matrix column index to reach the upper leftmost element

void FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::ReorthogonalizeVectors (ComplexVector* vectors, int nbrVectors, HermitianMatrix& matrix,
											      int rowShift, int columnShift)
{
  double TmpNorm = vectors[0].Norm();
  matrix(rowShift, columnShift) = TmpNorm;
  vectors[0] /= TmpNorm;
  for (int i = 1; i < nbrVectors; ++i)
    {
      MultipleComplexScalarProductOperation Operation (&(vectors[i]), 
						       vectors,   
						       i, this->TemporaryCoefficients);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < i; ++j)
	{
	  this->TemporaryCoefficients[j] = Conj(this->TemporaryCoefficients[j]);		  
	  matrix.SetMatrixElement(rowShift + i, columnShift + j, this->TemporaryCoefficients[j]);
	  this->TemporaryCoefficients[j] *= -1.0;
	}
      AddComplexLinearCombinationOperation Operation2 (&(vectors[i]), vectors, 
						       i, this->TemporaryCoefficients);	  
      Operation2.ApplyOperation(this->Architecture);
      TmpNorm = vectors[i].Norm();      
      matrix.SetMatrixElement(rowShift + i, columnShift + i, TmpNorm);
      vectors[i] /= TmpNorm;      
    }
}

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->EigenvaluePrecision);
  WriteLittleEndian(File, this->NbrEigenvalue);
  this->ReducedMatrix.WriteMatrix(File);
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    {
      WriteLittleEndian(File, this->PreviousWantedEigenvalues[i]);
    }
  File.close();
  return true;
}

// read current Lanczos state from disk
//
// return value = true if no error occurs

bool FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->EigenvaluePrecision);
  ReadLittleEndian(File, this->NbrEigenvalue);
  this->ReducedMatrix.ReadMatrix(File);
  if (this->PreviousWantedEigenvalues != 0)
    delete[] this->PreviousWantedEigenvalues;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    {
      ReadLittleEndian(File, this->PreviousWantedEigenvalues[i]);
    }
  File.close();  
  this->Diagonalize();
  return true;
}

// read a group of Lanczos vectors from disk
// 
// vectorArray = array where the vectors will be stored
// vectorAbsoluteIndex = absolute index of the first vector that has to be read from disk
// vectorRelativeIndex = index of the first vector that has to be read from disk within vectorArray
// nbrVectors = number of vectors to read from disk
// return value = true if no error occured

bool FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::ReadLanczosVectors(ComplexVector* vectorArray, int vectorAbsoluteIndex, 
											 int vectorRelativeIndex, int nbrVectors)
{
  char* TmpVectorName = new char [256];  
  for (int i = 0; i < nbrVectors; ++i)
    {
      sprintf(TmpVectorName, "vector.%d", vectorAbsoluteIndex + i);
      if (vectorArray[vectorRelativeIndex + i].ReadVector(TmpVectorName) == false)
	{
	  return false;
	}
   }
  delete[] TmpVectorName;
  return true;
}

// write a group of Lanczos vectors to disk
// 
// vectorArray = array where the vectors are stored
// vectorAbsoluteIndex = absolute index of the first vector that has to be written to disk
// vectorRelativeIndex = index of the first vector that has to be written to disk within vectorArray
// nbrVectors = number of vectors to written to disk
// return value = true if no error occured

bool FullReorthogonalizedComplexBlockLanczosAlgorithmWithDiskStorage::WriteLanczosVectors(ComplexVector* vectorArray, int vectorAbsoluteIndex, 
											  int vectorRelativeIndex, int nbrVectors)
{
  char* TmpVectorName = new char [256];  
  for (int i = 0; i < nbrVectors; ++i)
    {
      sprintf(TmpVectorName, "vector.%d", vectorAbsoluteIndex + i);
      if (vectorArray[vectorRelativeIndex + i].WriteVector(TmpVectorName) == false)
	{
	  return false;
	}
   }
  delete[] TmpVectorName;
  return true;
}
