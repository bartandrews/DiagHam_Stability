////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                 for non symmetric matrices using disk storage              //
//                                                                            //
//                        last modification : 07/01/2013                      //
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


#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
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
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
// leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// nbrTemporaryVectors = number of temporary that can be stored in memory
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

BasicArnoldiAlgorithmWithDiskStorage::BasicArnoldiAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter, 
									   bool highEnergy, bool leftFlag, bool resumeDiskFlag, int nbrTemporaryVectors, bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->NbrTemporaryVectors = nbrTemporaryVectors;
  this->ArnoldiVectors = new RealVector [this->NbrTemporaryVectors + 3];
  this->TemporaryCoefficients = new double [this->MaximumNumberIteration];
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = RealUpperHessenbergMatrix (this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix();
      this->ReducedMatrix = RealUpperHessenbergMatrix();
   }
  this->ResumeDiskFlag = resumeDiskFlag;
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

BasicArnoldiAlgorithmWithDiskStorage::BasicArnoldiAlgorithmWithDiskStorage(const BasicArnoldiAlgorithmWithDiskStorage& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->NbrTemporaryVectors = algorithm.NbrTemporaryVectors;
  this->ArnoldiVectors = new RealVector [3 + this->NbrTemporaryVectors];
  this->ResumeDiskFlag = algorithm.ResumeDiskFlag;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
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

BasicArnoldiAlgorithmWithDiskStorage::~BasicArnoldiAlgorithmWithDiskStorage() 
{
}

// initialize Lanczos algorithm with a random vector
//

void BasicArnoldiAlgorithmWithDiskStorage::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectors[0] = RealVector (Dimension);
  this->ArnoldiVectors[1] = RealVector (Dimension);
  this->ArnoldiVectors[2] = RealVector (Dimension);
  if (this->ResumeDiskFlag == false)
    {
      for (int i = 0; i < Dimension; i++)
	{
	  this->ArnoldiVectors[0][i] = (rand() - 32767) * 0.5;
	}
      this->ArnoldiVectors[0] /= this->ArnoldiVectors[0].Norm();
      this->Index = 0;
      this->TridiagonalizedMatrix.Resize(0, 0);
      this->ArnoldiVectors[0].WriteVector("vector.0"); 
    }
  else
    {
      this->ReadState();
    }
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicArnoldiAlgorithmWithDiskStorage::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if (this->ResumeDiskFlag == false)
    {
      this->ArnoldiVectors[0] = vector;
      this->ArnoldiVectors[1] = RealVector (Dimension);
      this->ArnoldiVectors[2] = RealVector (Dimension);
      this->Index = 0;
      this->TridiagonalizedMatrix.Resize(0, 0);
      this->ArnoldiVectors[0].WriteVector("vector.0"); 
    }
  else
    {
      this->ReadState();
    }
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& BasicArnoldiAlgorithmWithDiskStorage::GetGroundState()
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
  cout << "lapack is required for BasicArnoldiAlgorithmWithDiskStorage" << endl;
#endif
  if (this->HighEnergyFlag == false)
    SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector, true);
  else
    SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector, true);

  Complex* TmpCoefficents = new Complex [SortedDiagonalizedMatrix.GetNbrColumn()];
  for (int j = 0; j < SortedDiagonalizedMatrix.GetNbrColumn(); ++j)
    TmpCoefficents[j] = TmpEigenvector[0][j];
  int TmpTotalNbrReadVectors = 0;
  int TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, SortedDiagonalizedMatrix.GetNbrColumn());
  while (TmpNbrReadVectors > 0)
    {
      AddComplexLinearCombinationOperation Operation (&(this->GroundState), &(this->ArnoldiVectors[3]), 
						      TmpNbrReadVectors,  &(TmpCoefficents[TmpTotalNbrReadVectors]));
      Operation.ApplyOperation(this->Architecture);
      TmpTotalNbrReadVectors += TmpNbrReadVectors;
      TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, SortedDiagonalizedMatrix.GetNbrColumn());
    }
  delete[] TmpCoefficents;
  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* BasicArnoldiAlgorithmWithDiskStorage::GetEigenstates(int nbrEigenstates)
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
    {
      SortedDiagonalizedMatrix[i] = TmpDiag[i];
    }
#else
  cout << "lapack is required for BasicArnoldiAlgorithmWithDiskStorage" << endl;
#endif
  if (this->HighEnergyFlag == false)
    SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector, true);
  else
    SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector, true);

  Complex* TmpCoefficents = new Complex [SortedDiagonalizedMatrix.GetNbrColumn()];
  for (int i = 0; i < nbrEigenstates; ++i)
    {
      double TmpNorm = 1.0 / TmpEigenvector[i].Norm();
      for (int j = 0; j < SortedDiagonalizedMatrix.GetNbrColumn(); ++j)
	TmpCoefficents[j] = TmpEigenvector[i][j] * TmpNorm;
      Eigenstates[i] = ComplexVector(this->ArnoldiVectors[0].GetVectorDimension(), true);
      int TmpTotalNbrReadVectors = 0;
      int TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, SortedDiagonalizedMatrix.GetNbrColumn());
      while (TmpNbrReadVectors > 0)
	{
	  AddComplexLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->ArnoldiVectors[3]), 
							  TmpNbrReadVectors,  &(TmpCoefficents[TmpTotalNbrReadVectors]));
	  Operation.ApplyOperation(this->Architecture);
	  TmpTotalNbrReadVectors += TmpNbrReadVectors;
	  TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, SortedDiagonalizedMatrix.GetNbrColumn());
	}
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicArnoldiAlgorithmWithDiskStorage::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      if (nbrIter < 2)
	nbrIter = 2;
      Dimension = nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);

      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->ArnoldiVectors[0]), &(this->ArnoldiVectors[1]));
      Operation1.ApplyOperation(this->Architecture);
      this->ReducedMatrix.SetMatrixElement(0, 0, (this->ArnoldiVectors[0] * this->ArnoldiVectors[1]));
      double Tmp;
      this->ReducedMatrix.GetMatrixElement(0, 0, Tmp);
      this->ArnoldiVectors[1].AddLinearCombination(-Tmp, this->ArnoldiVectors[0]);
      this->ReducedMatrix.SetMatrixElement(1, 0, this->ArnoldiVectors[1].Norm()); 
      this->ReducedMatrix.GetMatrixElement(1, 0, Tmp);
      this->ArnoldiVectors[1] /=  Tmp; 
      this->ArnoldiVectors[1].WriteVector("vector.1");
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->ArnoldiVectors[1]), &(this->ArnoldiVectors[2]));
      Operation2.ApplyOperation(this->Architecture);
      this->ReducedMatrix.SetMatrixElement(0, 1, (this->ArnoldiVectors[0] * this->ArnoldiVectors[2]));
      this->ReducedMatrix.SetMatrixElement(1, 1, (this->ArnoldiVectors[1] * this->ArnoldiVectors[2]));
    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter;
      this->ReducedMatrix.ResizeAndClean(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; ++i)
    {
      if (this->ResumeDiskFlag == false)
	{
	  for (int k = 0; k < i; ++k)
	    {
	      this->ReducedMatrix.GetMatrixElement(k, i - 1, this->TemporaryCoefficients[k]);
	      this->TemporaryCoefficients[k] *= -1.0;
	    }
	  int TmpTotalNbrReadVectors = 0;
	  int TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, i);
	  while (TmpNbrReadVectors > 0)
	    {
	      AddRealLinearCombinationOperation Operation (&(this->ArnoldiVectors[2]), &(this->ArnoldiVectors[3]), 
							   TmpNbrReadVectors, &(this->TemporaryCoefficients[TmpTotalNbrReadVectors]));
	      Operation.ApplyOperation(this->Architecture);
	      TmpTotalNbrReadVectors += TmpNbrReadVectors;
	      TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, i);
	    }
	  double VectorNorm = this->ArnoldiVectors[2].Norm();
	  this->ReducedMatrix.SetMatrixElement(i, i - 1, VectorNorm);
	  if (VectorNorm < 1e-5)
	    {
	      cout << "subspace !!! " << i << endl;
	    }
	  this->ArnoldiVectors[2] /= VectorNorm;
	  char* TmpVectorName = new char [256];
	  sprintf(TmpVectorName, "vector.%d", i);
	  this->ArnoldiVectors[2].WriteVector(TmpVectorName);
	  delete[] TmpVectorName;
	  this->WriteState();
	}
      else
	{
	  this->ResumeDiskFlag = false;
	}

      RealVector TmpVector = this->ArnoldiVectors[0];
      this->ArnoldiVectors[0] = this->ArnoldiVectors[1];
      this->ArnoldiVectors[1] = this->ArnoldiVectors[2];
      this->ArnoldiVectors[2] = TmpVector;

      this->Index++;
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->ArnoldiVectors[1]), &(this->ArnoldiVectors[2]));
      Operation.ApplyOperation(this->Architecture);
      int TmpTotalNbrReadVectors = 0;
      int TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, i + 1);
      while (TmpNbrReadVectors > 0)
	{
	  MultipleRealScalarProductOperation Operation2 (&(this->ArnoldiVectors[2]), &(this->ArnoldiVectors[3]), 
							 TmpNbrReadVectors,  &(this->TemporaryCoefficients[TmpTotalNbrReadVectors]));
	  Operation2.ApplyOperation(this->Architecture);
	  TmpTotalNbrReadVectors += TmpNbrReadVectors;
	  TmpNbrReadVectors = this->ReadTemporaryVectors(TmpTotalNbrReadVectors, i + 1);
	}
      for (int j = 0; j <= i; ++j)
	{
	  this->ReducedMatrix.SetMatrixElement(j, i, this->TemporaryCoefficients[j]);
	}
    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->ComplexPreviousWantedEigenvalues[i] = this->ComplexDiagonalizedMatrix[i];
      this->Diagonalize();
      if (this->HighEnergyFlag == false)
	this->ComplexDiagonalizedMatrix.SortMatrixUpOrder(true);
      else
	this->ComplexDiagonalizedMatrix.SortMatrixDownOrder(true);
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

  
// write current Lanczos state on disk
//
// return value = true if no error occurs

bool BasicArnoldiAlgorithmWithDiskStorage::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->NbrEigenvalue);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue.Re);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue.Im);
  for (int i = 0; i < this->NbrEigenvalue; ++i)    
    {    
      WriteLittleEndian(File, this->ComplexPreviousWantedEigenvalues[i].Re);
      WriteLittleEndian(File, this->ComplexPreviousWantedEigenvalues[i].Im);
    }
  this->ReducedMatrix.WriteMatrix(File);
  File.close();  
  return true;
}

// read current Lanczos state from disk
//
// return value = true if no error occurs

bool BasicArnoldiAlgorithmWithDiskStorage::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->NbrEigenvalue);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue.Re);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue.Im);
  for (int i = 0; i < this->NbrEigenvalue; ++i)    
    {    
      ReadLittleEndian(File, this->ComplexPreviousWantedEigenvalues[i].Re);
      ReadLittleEndian(File, this->ComplexPreviousWantedEigenvalues[i].Im);
    }
  this->ReducedMatrix.ReadMatrix(File);
  File.close();  
  char* TmpVectorName = new char [256];
  sprintf(TmpVectorName, "vector.%d", this->Index);
  this->ArnoldiVectors[0].ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 1));
  this->ArnoldiVectors[1].ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 2));
  this->ArnoldiVectors[2].ReadVector(TmpVectorName);
  delete[] TmpVectorName;
  return true;
}

// read several temporary vectors stored o disk
//
// firstVector = index of the first vector to read
// totalNbrVectors = total number of temporary vectors
// return value = number of vectors that have been read

int BasicArnoldiAlgorithmWithDiskStorage::ReadTemporaryVectors(int firstVector, int totalNbrVectors)
{
  char* TmpVectorName = new char [256];
  int TmpNbrVector = (totalNbrVectors - firstVector);
  if (TmpNbrVector > this->NbrTemporaryVectors)
    TmpNbrVector = this->NbrTemporaryVectors;
  for (int i = 0; i < TmpNbrVector; ++i)
    {
      sprintf(TmpVectorName, "vector.%d", (firstVector + i));
      this->ArnoldiVectors[3 + i].ReadVector(TmpVectorName);
    }
  delete[] TmpVectorName;
  return TmpNbrVector;
}
