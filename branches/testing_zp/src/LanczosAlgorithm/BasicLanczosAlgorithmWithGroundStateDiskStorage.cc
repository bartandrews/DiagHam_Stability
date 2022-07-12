////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2006 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of basic Lanczos algorithm with real vectors           //
//            and ground state evaluation and disk storage capabilities       //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 10/01/2006                      //
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


#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundStateDiskStorage.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>


using std::cout;
using std::endl;
using std::ios;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrIter = maximum number of iterations when evaluating the ground state eigenvector (0 if all iterations needed for convergence have to be done)
// maxIter = an approximation of maximal number of iteration

BasicLanczosAlgorithmWithGroundStateDiskStorage::BasicLanczosAlgorithmWithGroundStateDiskStorage(AbstractArchitecture* architecture,
												 int nbrIter, int maxIter) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = RealVector();
  this->V2 = RealVector();
  this->V3 = RealVector();
  this->InitialState = RealVector();
  this->GroundState = RealVector();
  this->GroundStateFlag = false;
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
  this->NbrIterationsGroundState = nbrIter;
  this->Architecture = architecture;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->NbrEigenvalue = 1;
  this->GroundStateEvaluationFlag = 0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

BasicLanczosAlgorithmWithGroundStateDiskStorage::BasicLanczosAlgorithmWithGroundStateDiskStorage(const BasicLanczosAlgorithmWithGroundStateDiskStorage& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->V3 = algorithm.V3;
  this->InitialState = algorithm.InitialState;
  this->GroundState = algorithm.GroundState;
  this->GroundStateFlag = algorithm.GroundStateFlag;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->NbrEigenvalue = 1;
  this->GroundStateEvaluationFlag = algorithm.GroundStateEvaluationFlag;
  this->NbrIterationsGroundState = algorithm.NbrIterationsGroundState;
}

// destructor
//

BasicLanczosAlgorithmWithGroundStateDiskStorage::~BasicLanczosAlgorithmWithGroundStateDiskStorage() 
{
}

// initialize Lanczos algorithm with a random vector
//

void BasicLanczosAlgorithmWithGroundStateDiskStorage::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = RealVector (Dimension);
  this->V2 = RealVector (Dimension);
  this->V3 = RealVector (Dimension);
  int Shift = RAND_MAX / 2;
  double Scale = 1.0 / ((double) Shift);
  for (int i = 0; i < Dimension; i++)
    {
      this->V1[i] = Scale * ((double) (rand() - Shift));
    }
  this->V1 /= this->V1.Norm();
  this->InitialState = RealVector (this->V1, true);
  this->InitialState.WriteVector("vector.0");
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicLanczosAlgorithmWithGroundStateDiskStorage::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = vector;
  this->V2 = RealVector (Dimension);
  this->V3 = RealVector (Dimension);
  this->InitialState = RealVector (vector, true);
  this->Index = 0;
  this->InitialState.WriteVector("vector.0");
  this->GroundStateFlag = false;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// resume Lanczos algorithm from disk datas in current directory
//

void BasicLanczosAlgorithmWithGroundStateDiskStorage::ResumeLanczosAlgorithm()
{ 
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ReadState();
  this->InitialState = RealVector (Dimension);
  this->V1 = RealVector (Dimension);
  this->V2 = RealVector (Dimension);
  this->V3 = RealVector (Dimension);
  this->InitialState.ReadVector("vector.0");
  this->V1.ReadVector("vector.1");
  this->V2.ReadVector("vector.2");
  this->V3.ReadVector("vector.3");
  if (this->TestConvergence() == true)
    {
      this->GroundState.ReadVector("vector.4");
    }
}
  
// get the n first eigenstates (limited to the ground state for this class, return NULL if nbrEigenstates > 1)
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* BasicLanczosAlgorithmWithGroundStateDiskStorage::GetEigenstates(int nbrEigenstates)
{
  if (nbrEigenstates != 1)
    {
      return 0;
    }
  else
    {
      this->GetGroundState();
      RealVector* TmpVectors = new RealVector [1];
      TmpVectors[0] = this->GroundState;
      return TmpVectors;
    }
}

// get last produced vector
//
// return value = reference on last produced vector

Vector& BasicLanczosAlgorithmWithGroundStateDiskStorage::GetGroundState()
{
  if (this->GroundStateFlag == false)
    {
      RealMatrix TmpEigenvector (this->TridiagonalizedMatrix.GetNbrRow(), this->TridiagonalizedMatrix.GetNbrRow(), true);
      for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
	TmpEigenvector(i, i) = 1.0;
      
      RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->TridiagonalizedMatrix.GetNbrRow());
      SortedDiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
      SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
      SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
      double* TmpComponents = new double [this->TridiagonalizedMatrix.GetNbrRow()];
      for (int j = 0; j < this->TridiagonalizedMatrix.GetNbrRow(); ++j)
	{
	  TmpComponents[j] = TmpEigenvector(j, 0);
	}

      if (this->GroundStateEvaluationFlag == 0)
	{
	  this->GroundState.Copy(this->InitialState, TmpComponents[0]);
	  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->InitialState, &this->V3);
	  Operation1.ApplyOperation(this->Architecture);
	  this->V3.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(0), this->InitialState);
	  this->V3 /= this->V3.Norm();
	  this->V2.Copy(this->InitialState);
	  this->GroundState.AddLinearCombination(TmpComponents[1], this->V3);
	  if (this->Index > 0)
	    {
	      this->NbrIterationsGroundState -= (this->Index % this->NbrIterationsGroundState);
	      if (this->NbrIterationsGroundState <= 0)
		{
		  this->NbrIterationsGroundState = 1;
		}
	    }
	  this->Index = 2;
	  this->GroundStateEvaluationFlag = 1;
	  this->WriteState();
	}
      int Lim = this->DiagonalizedMatrix.GetNbrRow();
      if (this->NbrIterationsGroundState > 0)
	{
	  cout << this->NbrIterationsGroundState << " " << this->Index << " " << Lim << " " << this->DiagonalizedMatrix.GetNbrRow() << endl;
	  Lim = this->Index + this->NbrIterationsGroundState;
	  if (Lim > this->DiagonalizedMatrix.GetNbrRow())
	    Lim = this->DiagonalizedMatrix.GetNbrRow();
	}
      for (int i = this->Index; i < Lim; ++i)
	{
	  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V3, &this->V1);
	  Operation1.ApplyOperation(this->Architecture);
	  this->V1.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(i - 1), this->V3, 
					-this->TridiagonalizedMatrix.UpperDiagonalElement(i - 2), this->V2);
	  this->V1 /= this->V1.Norm();
	  this->GroundState.AddLinearCombination(TmpComponents[i], this->V1);
	  RealVector TmpV (this->V2);
	  this->V2 = this->V3;
	  this->V3 = this->V1;
	  this->V1 = TmpV;
	  this->V1.WriteVector("vector.1");
	  this->V2.WriteVector("vector.2");
	  this->V3.WriteVector("vector.3");
	  this->GroundState.WriteVector("vector.4");
	  cout << i << "/" << Lim << "           \r";
	  cout.flush();
	  ++this->Index;
	  this->WriteState();
	}
      cout << endl;
      if (Lim != this->DiagonalizedMatrix.GetNbrRow())
	{
	  cout << "need " << (this->DiagonalizedMatrix.GetNbrRow() - Lim) << " more iterations to get the ground state" << endl;
	}
      this->GroundState /= this->GroundState.Norm();
      this->GroundStateFlag = true;
      delete[] TmpComponents;
    }
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicLanczosAlgorithmWithGroundStateDiskStorage::RunLanczosAlgorithm (int nbrIter) 
{
  this->GroundStateFlag = false;
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V1, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->V1 * this->V2);
      this->V2.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index), 
				    this->V1);
      this->V2 /= this->V2.Norm(); 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &this->V2, &this->V3);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3);
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  RealVector* TmpVector = new RealVector[2];
  double* TmpCoefficient = new double[2];
  RealVector* TmpVectorScalarProduct[2];
  double TmpScalarProduct[2];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      TmpVector[0] = this->V1;
      TmpVector[1] = this->V2;
      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
      AddRealLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);
      this->V3 /= this->V3.Norm();
      RealVector TmpV (this->V1);
      this->V1 = this->V2;
      this->V2 = this->V3;
      this->V3 = TmpV;
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V2, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
      TmpVectorScalarProduct[0] = &(this->V1);
      TmpVectorScalarProduct[1] = &(this->V2);
      MultipleRealScalarProductOperation Operation2 (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = TmpScalarProduct[0];
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = TmpScalarProduct[1];
      this->V1.WriteVector("vector.1");
      this->V2.WriteVector("vector.2");
      this->V3.WriteVector("vector.3");
      this->WriteState();
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
  this->WriteState();
}
  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool BasicLanczosAlgorithmWithGroundStateDiskStorage::TestConvergence ()
{
  if ((fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
       (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1)))))
    return true;
  else
    return false;
}

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool BasicLanczosAlgorithmWithGroundStateDiskStorage::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->GroundStateEvaluationFlag);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->EigenvaluePrecision);
  WriteLittleEndian(File, this->NbrEigenvalue);
  int TmpDimension = this->TridiagonalizedMatrix.GetNbrRow();
  WriteLittleEndian(File, TmpDimension);
  for (int i = 0; i < TmpDimension; ++i)    
    {    
      WriteLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  --TmpDimension;
  for (int i = 0; i < TmpDimension; ++i)
    {
      WriteLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  File.close();
  return true;
}

// read current Lanczos state from disk
//
// return value = true if no error occurs

bool BasicLanczosAlgorithmWithGroundStateDiskStorage::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->GroundStateEvaluationFlag);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->EigenvaluePrecision);
  ReadLittleEndian(File, this->NbrEigenvalue);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
  for (int i = 0; i < TmpDimension; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  --TmpDimension;
  for (int i = 0; i < TmpDimension; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  File.close();  
  this->Diagonalize();
  this->DiagonalizedMatrix.SortMatrixUpOrder();
  return true;
}

