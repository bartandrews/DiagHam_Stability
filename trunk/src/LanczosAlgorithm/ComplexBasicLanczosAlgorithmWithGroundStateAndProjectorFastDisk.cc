////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of basic Lanczos algorithm with real vectors           //
//                         and ground state evaluation                        //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 17/09/2002                      //
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


#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;



// default constructor
//
// nbrProjectors = dimension of the projector subspace
// projectorVectors = array that contains the vectors that spans the projector subspace
// projectorCoefficient = energy scale in front of the projector
// indexShiftFlag = true if the eigenstate indices have to be shifted
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// returnAtConvergenceFlag = flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrProjectors, ComplexVector* projectorVectors, 
																 double projectorCoefficient, bool indexShiftFlag, 
																 AbstractArchitecture* architecture, 
																 int maxIter, bool diskFlag, bool resumeDiskFlag, bool returnAtConvergenceFlag) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->V3 = ComplexVector();
  this->InitialState = ComplexVector();
  this->GroundState = ComplexVector();
  this->GroundStateFlag = false;
  this->DiskFlag = diskFlag;
  this->ResumeDiskFlag = resumeDiskFlag;
  this->NbrProjectors = nbrProjectors;
  this->InitialNbrProjectors = nbrProjectors;
  this->ProjectorVectors = new ComplexVector [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    this->ProjectorVectors[i] = projectorVectors[i];
  this->ProjectorCoefficient = projectorCoefficient;
  this->IndexShiftFlag = indexShiftFlag;
  this->AutomaticProjectorConstructionFlag = false;
  this->ProjectorEigenvalues = 0;
  this->ReturnAtConvergenceFlag = returnAtConvergenceFlag;
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
  this->NbrEigenvalue = 1;
}

// constructor using automatic projector construction 
//
// nbrEigenvalues = number of eigenvalues/eigenstates to compute
// projectorCoefficient = energy scale in front of the projector
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// returnAtConvergenceFlag = flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrEigenvalues, double projectorCoefficient,
																 AbstractArchitecture* architecture, 
																 int maxIter, bool diskFlag, bool resumeDiskFlag, bool returnAtConvergenceFlag)
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->V3 = ComplexVector();
  this->InitialState = ComplexVector();
  this->GroundState = ComplexVector();
  this->GroundStateFlag = false;
  this->DiskFlag = diskFlag;
  this->ResumeDiskFlag = resumeDiskFlag;
  this->NbrEigenvalue = nbrEigenvalues;
  this->NbrProjectors = 0;
  this->InitialNbrProjectors = 0;
  this->ProjectorVectors = new ComplexVector [this->NbrEigenvalue - 1];
  this->ProjectorEigenvalues = new double [this->NbrEigenvalue - 1]; 
  this->ProjectorCoefficient = projectorCoefficient;
  this->IndexShiftFlag = false;
  this->AutomaticProjectorConstructionFlag = true;
  this->ReturnAtConvergenceFlag = returnAtConvergenceFlag;
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

// constructor using both automatic projector construction and an initial set of projectors
//
// nbrEigenvalues = number of eigenvalues/eigenstates to compute
// nbrProjectors = dimension of the projector subspace
// projectorVectors = array that contains the vectors that spans the projector subspace
// projectorCoefficient = energy scale in front of the projector
// indexShiftFlag = true if the eigenstate indices have to be shifted
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// returnAtConvergenceFlag = flag indicating whether the algorithm should return before replaying the Lanczos run (useful for large parallel runs)
ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(int nbrEigenvalues,
																 int nbrProjectors, ComplexVector* projectorVectors,
																 double projectorCoefficient, bool indexShiftFlag,
																 AbstractArchitecture* architecture, 
																 int maxIter, bool diskFlag, bool resumeDiskFlag, bool returnAtConvergenceFlag)
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->V3 = ComplexVector();
  this->InitialState = ComplexVector();
  this->GroundState = ComplexVector();
  this->GroundStateFlag = false;
  this->DiskFlag = diskFlag;
  this->ResumeDiskFlag = resumeDiskFlag;
  this->NbrEigenvalue = nbrEigenvalues;
  this->NbrProjectors = nbrProjectors;
  this->InitialNbrProjectors = nbrProjectors;
  this->ProjectorVectors = new ComplexVector [this->NbrEigenvalue - 1 + this->InitialNbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    this->ProjectorVectors[i] = projectorVectors[i];
  this->ProjectorEigenvalues = new double [this->NbrEigenvalue - 1]; 
  this->ProjectorCoefficient = projectorCoefficient;
  this->IndexShiftFlag = indexShiftFlag;
  this->AutomaticProjectorConstructionFlag = true;
  this->ReturnAtConvergenceFlag = returnAtConvergenceFlag;
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

ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk(const ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->V3 = algorithm.V3;
  this->DiskFlag = algorithm.DiskFlag;
  this->ResumeDiskFlag = algorithm.ResumeDiskFlag;
  this->InitialState = algorithm.InitialState;
  this->GroundState = algorithm.GroundState;
  this->GroundStateFlag = algorithm.GroundStateFlag;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->AutomaticProjectorConstructionFlag = algorithm.AutomaticProjectorConstructionFlag;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->NbrProjectors = algorithm.NbrProjectors;
  this->InitialNbrProjectors = algorithm.InitialNbrProjectors;
  this->ReturnAtConvergenceFlag = algorithm.ReturnAtConvergenceFlag;
  if (this->AutomaticProjectorConstructionFlag == false)
    {
      this->ProjectorVectors = new ComplexVector [this->NbrProjectors];
      for (int i = 0; i < this->NbrProjectors; ++i)
	this->ProjectorVectors[i] = algorithm.ProjectorVectors[i];
      this->ProjectorEigenvalues = 0;
    }
  else
    {
      this->ProjectorVectors = new ComplexVector [this->NbrEigenvalue - 1 + this->InitialNbrProjectors];     
      for (int i = 0; i < this->NbrProjectors; ++i)
	this->ProjectorVectors[i] = algorithm.ProjectorVectors[i];
      this->ProjectorEigenvalues = new double [this->NbrEigenvalue - 1]; 
      for (int i = 0; i < this->NbrProjectors; ++i)
	this->ProjectorEigenvalues[i] = algorithm.ProjectorEigenvalues[i];
    }
  this->ProjectorCoefficient = algorithm.ProjectorCoefficient;
  this->IndexShiftFlag = algorithm.IndexShiftFlag;
}

// destructor
//

ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::~ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk() 
{
  delete[] this->ProjectorVectors;
  if (this->ProjectorEigenvalues != 0)
    delete[] this->ProjectorEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if ((this->AutomaticProjectorConstructionFlag == false) || ((this->NbrProjectors - this->InitialNbrProjectors) == 0) || (this->ResumeDiskFlag == true))
    {
      this->V1 = ComplexVector (Dimension);
      this->V2 = ComplexVector (Dimension);
      this->V3 = ComplexVector (Dimension);
    }
  if (this->ResumeDiskFlag == false)
    {
      int Shift = RAND_MAX / 2;
      double Scale = 1.0 / ((double) Shift);
      for (int i = 0; i < Dimension; i++)
	{
	  this->V1.Re(i) = Scale * ((double) (rand() - Shift));
	  this->V1.Im(i) = Scale * ((double) (rand() - Shift));
	}
      this->V1 /= this->V1.Norm();
      if (this->NbrProjectors > 0)
	{
	  Complex* TmpScalarProduct = new Complex[this->NbrProjectors];
	  MultipleComplexScalarProductOperation Operation1 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
	  Operation1.ApplyOperation(this->Architecture);	
	  for (int i = 0; i < this->NbrProjectors; ++i)
	    TmpScalarProduct[i] = -Conj(TmpScalarProduct[i]);
	  AddComplexLinearCombinationOperation Operation2 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
	  Operation2.ApplyOperation(this->Architecture);
	  delete[] TmpScalarProduct;
	  this->V1 /= this->V1.Norm();
	}
      if (this->DiskFlag == false)
	this->InitialState = ComplexVector (this->V1, true);
      else
	this->V1.WriteVector("vector.0");
      this->Index = 0;
      this->TridiagonalizedMatrix.Resize(0, 0);
    }
  else
    {
      this->ReadState();
    }
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if (this->ResumeDiskFlag == false)
    {
      this->V1 = ((ComplexVector &) vector);
      this->V2 = ComplexVector (Dimension);
      this->V3 = ComplexVector (Dimension);
      if (this->NbrProjectors > 0)
	{
	  Complex* TmpScalarProduct = new Complex[this->NbrProjectors];
	  MultipleComplexScalarProductOperation Operation1 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
	  Operation1.ApplyOperation(this->Architecture);	
	  for (int i = 0; i < this->NbrProjectors; ++i)
	    TmpScalarProduct[i] = -Conj(TmpScalarProduct[i]);
	  AddComplexLinearCombinationOperation Operation2 (&(this->V1), this->ProjectorVectors, this->NbrProjectors, TmpScalarProduct);
	  Operation2.ApplyOperation(this->Architecture);
	  delete[] TmpScalarProduct;
	  this->V1 /= this->V1.Norm();
	}
      if (this->DiskFlag == false)
	this->InitialState = ComplexVector (vector);
      else
	this->V1.WriteVector("vector.0");
      this->Index = 0;
      this->GroundStateFlag = false;
      this->TridiagonalizedMatrix.Resize(0, 0);
    }
  else
    {
      this->V1 = ComplexVector (Dimension);
      this->V2 = ComplexVector (Dimension);
      this->V3 = ComplexVector (Dimension);
      this->ReadState();
    }
}

// get the n first eigenstates (limited to the ground state from this class, return NULL if nbrEigenstates > 1)
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::GetEigenstates(int nbrEigenstates)
{
  if (this->AutomaticProjectorConstructionFlag == false)
    {
      if (nbrEigenstates != 1)
	{
	  return 0;
	}
      else
	{
	  this->GetGroundState();
	  ComplexVector* TmpVectors = new ComplexVector[1];
	  TmpVectors[0] = this->GroundState;
	  return TmpVectors;
	}
    }
  else
    {
      if (nbrEigenstates > this->NbrEigenvalue)
	{
	  return 0;	  
	}
      ComplexVector* TmpVectors = new ComplexVector[nbrEigenstates];
      if (nbrEigenstates < this->NbrEigenvalue)
	{
	  for (int i = 0; i < nbrEigenstates; ++i)
	    TmpVectors[i] = this->ProjectorVectors[this->InitialNbrProjectors + i];
	}
      else
	{
	  for (int i = this->InitialNbrProjectors; i < this->NbrProjectors; ++i)
	    TmpVectors[i - this->InitialNbrProjectors] = this->ProjectorVectors[i];
	  this->GetGroundState();
	  TmpVectors[this->NbrEigenvalue - 1] = this->GroundState;
	}
      this->NbrProjectors = this->InitialNbrProjectors;
      return TmpVectors;    
    } 
}

// get last produced vector
//
// return value = reference on last produced vector

Vector& ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::GetGroundState()
{
  if (this->GroundStateFlag == false)
    {
      // changed from ComplexMatrix
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

      if (this->DiskFlag == false)
	{
	  // save initial vector and current state for possible recovery.
	  this->InitialState.WriteVector("initialvector");
	  this->WriteState();
	  //
	  double* TmpCoefficient = new double[2];
	  this->GroundState.Copy(this->InitialState, TmpComponents[0]);
	  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->InitialState, &this->V3);
	  Operation1.ApplyOperation(this->Architecture);
	  this->AddProjectorContribution(this->InitialState, this->V3);
	  this->V3.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(0), this->InitialState);
	  this->V3 /= this->V3.Norm();
	  this->V2.Copy(this->InitialState);
	  this->GroundState.AddLinearCombination(TmpComponents[1], this->V3);
	  for (int i = 2; i < this->DiagonalizedMatrix.GetNbrRow(); ++i)
	    {
	      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V3, &this->V1);
	      Operation1.ApplyOperation(this->Architecture);
	      this->AddProjectorContribution(this->V3, this->V1);
  	      ComplexVector* TmpVector = new ComplexVector[2];
	      TmpVector[0] = this->V2;
	      TmpVector[1] = this->V3;
	      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(i - 2);
	      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(i - 1);
	      AddComplexLinearCombinationOperation Operation4 (&(this->V1),  TmpVector, 2, TmpCoefficient);
	      Operation4.ApplyOperation(this->Architecture);
	      delete[] TmpVector;
	      this->V1 /= this->V1.Norm();
	      this->GroundState.AddLinearCombination(TmpComponents[i], this->V1);
	      ComplexVector TmpV (this->V2);
	      this->V2 = this->V3;
	      this->V3 = this->V1;
	      this->V1 = TmpV;
	      cout << i << "/" << this->DiagonalizedMatrix.GetNbrRow() << "           \r";
 	      cout.flush();
	    }
	}
      else
	{ 
	  if (this->ReturnAtConvergenceFlag)
	    {
	      cout << "returning prior to calculation of eigenstate: state can be obtained with ReplayFastLanczos."<<endl;
	      exit(0);
	    }
	  this->V1.ReadVector("vector.0");	      
	  this->GroundState.Copy(this->V1, TmpComponents[0]);
	  char* TmpVectorName = new char [256];
	  for (int i = 1; i < this->DiagonalizedMatrix.GetNbrRow(); ++i)
	    {
	      sprintf(TmpVectorName, "vector.%d", i);
	      this->V1.ReadVector(TmpVectorName);	      
	      this->GroundState.AddLinearCombination(TmpComponents[i], this->V1);	      
	      cout << i << "/" << this->DiagonalizedMatrix.GetNbrRow() << "           \r";
	      cout.flush();
	    }	  
	  delete[] TmpVectorName;
	}
      cout << endl;
      this->GroundState /= this->GroundState.Norm();
      this->GroundStateFlag = true;
      delete[] TmpComponents;
    }
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::RunLanczosAlgorithm (int nbrIter) 
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
      this->AddProjectorContribution(this->V1, this->V2);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->V1 * this->V2).Re;
      this->V2.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index), this->V1);
      this->V2 /= this->V2.Norm(); 
      if (this->DiskFlag == true)
	this->V2.WriteVector("vector.1");
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
  Complex* TmpCoefficient = new Complex[2];
  Complex TmpScalarProduct[2];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      ComplexVector* TmpVector = new ComplexVector[2];
      if (this->ResumeDiskFlag == false)
	{
	  TmpVector[0] = this->V1;
	  TmpVector[1] = this->V2;
	  TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
	  TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
	  AddComplexLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
	  Operation4.ApplyOperation(this->Architecture);
	  delete[] TmpVector;
	  this->V3 /= this->V3.Norm();
	  if (this->DiskFlag == true)
	    {
	      char* TmpVectorName = new char [256];
	      sprintf(TmpVectorName, "vector.%d", i);
	      this->V3.WriteVector(TmpVectorName);
	      delete[] TmpVectorName;
	      this->WriteState();
	    }
	}
      else
	{
	  this->ResumeDiskFlag = false;
	}
      if (this->DiskFlag == true)
	{
	  ComplexVector TmpV (this->V2);
	  this->V2 = this->V3;
	  this->V3 = TmpV;	  
	  this->V1 = ComplexVector();
	}
      else
	{
	  ComplexVector TmpV (this->V1);
	  this->V1 = this->V2;
	  this->V2 = this->V3;
	  this->V3 = TmpV;
	}
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V2, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
      this->AddProjectorContribution(this->V2, this->V3);
      if (this->DiskFlag == true)
	{
	  char* TmpVectorName = new char [256];
	  sprintf(TmpVectorName, "vector.%d", (i - 1));
	  this->V1.ReadVector(TmpVectorName);
	  delete[] TmpVectorName;
	}      
      ComplexVector* TmpVectorScalarProduct[2];
      TmpVectorScalarProduct[0] = &(this->V1);
      TmpVectorScalarProduct[1] = &(this->V2);
      MultipleComplexScalarProductOperation Operation2 (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = TmpScalarProduct[0].Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = TmpScalarProduct[1].Re;
    }
  delete[] TmpCoefficient;
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = this->DiagonalizedMatrix.DiagonalElement(0);//this->NbrEigenvalue - 1);
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * this->DiagonalizedMatrix.DiagonalElement(0);//this->NbrEigenvalue - 1);
    }
}
  
// optional shift of the eigenstate file name indices
//
// return value = index shift

int ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::EigenstateIndexShift()
{
  if (this->IndexShiftFlag == true)
    {
      return this->InitialNbrProjectors;
    }
  else
    {
      return 0;
    }
}

// add the projector contribution to the hamiltonian-vector multiplication
//
// initialVector = reference on the initial vector
// destinationVector = reference on the destination vector 

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::AddProjectorContribution(ComplexVector& initialVector, ComplexVector& destinationVector)
{
  if (this->NbrProjectors > 0)
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
}

// test if convergence has been reached
//
// return value = true if convergence has been reached

bool ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::TestConvergence ()
{
  if ((fabs(this->DiagonalizedMatrix.DiagonalElement(0) - this->PreviousLastWantedEigenvalue) < 
       (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(0)))))
    {
      if (this->AutomaticProjectorConstructionFlag == true)
	{
	  if ((this->NbrProjectors - this->InitialNbrProjectors+ 1) == this->NbrEigenvalue)
	    return true;
	  else
	    {
	      cout << "found eigenstate " << (this->NbrProjectors - this->InitialNbrProjectors) <<endl;
	      cout << "computing eigenstate " << (this->NbrProjectors - this->InitialNbrProjectors) << endl;
	      this->GetGroundState();
	      this->ProjectorVectors[this->NbrProjectors] = ComplexVector(this->GroundState, true);
	      if (this->NbrProjectors > 0)
		{
		  cout << "checking orthogonality with previously other projector states" << endl;
		  for (int i = 0; i < this->NbrProjectors; ++i)
		    {
		      cout << "< " << i << " | " << this->NbrProjectors << " > = " 
			   << (this->ProjectorVectors[i] * this->ProjectorVectors[this->NbrProjectors]) << endl;
		    }
		}
	      this->ProjectorEigenvalues[this->NbrProjectors - this->InitialNbrProjectors] = this->DiagonalizedMatrix.DiagonalElement(0);
	      if (this->DiskFlag == true)
		{
		  char* TmpVectorName = new char [256];
		  sprintf(TmpVectorName, "vectorprojector.%d", this->NbrProjectors);
		  this->ProjectorVectors[this->NbrProjectors].WriteVector(TmpVectorName);
		  delete[] TmpVectorName;
		}
	      this->NbrProjectors++;
	      this->GroundStateFlag = false;
	      cout << "reinitializing Lanczos " << endl;
	      this->Index = 0;
	      this->PreviousLastWantedEigenvalue = 0.0;	  
	      
	      this->TridiagonalizedMatrix.ClearMatrix();
	      this->DiagonalizedMatrix.ClearMatrix();
	      this->InitializeLanczosAlgorithm();
	      cout << "starting Lanczos " << endl;
	      return false;
	    }
	}
      else
	{
	  return true;
	}
    }
  else
    return false;
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::GetEigenvalues (double*& eigenvalues, int nbrEigenvalues)
{
  if (this->AutomaticProjectorConstructionFlag == false)
    {
      eigenvalues = new double [nbrEigenvalues];
      for (int i = 0; i < nbrEigenvalues; ++i)
	{
	  eigenvalues[i] = this->DiagonalizedMatrix(i, i);
	}
    }
  else
    {
      eigenvalues = new double [nbrEigenvalues];
      for (int i = 0; i < (this->NbrProjectors - this->InitialNbrProjectors); ++i)
	{
	  eigenvalues[i] = this->ProjectorEigenvalues[i];
	}
      for (int i = this->NbrProjectors; i < nbrEigenvalues; ++i)
	{
	  eigenvalues[i] = this->DiagonalizedMatrix(i, i);
	}      
    }
}

// get current diagonalized matrix
//
// return value = reference on current diagonalized matrix

RealTriDiagonalSymmetricMatrix& ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::GetDiagonalizedMatrix ()
{
  if (this->AutomaticProjectorConstructionFlag == false)
    {
      return this->DiagonalizedMatrix;
    }
  else
    {
      this->FullDiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->DiagonalizedMatrix.GetNbrRow() + (this->NbrProjectors - this->InitialNbrProjectors), true);
      for (int i = 0; i < (this->NbrProjectors - this->InitialNbrProjectors); ++i)
	{
	  this->FullDiagonalizedMatrix(i, i) = this->ProjectorEigenvalues[i];
	}
      for (int i = 0; i < this->DiagonalizedMatrix.GetNbrRow(); ++i)
	{
	  this->FullDiagonalizedMatrix(i + this->NbrProjectors - this->InitialNbrProjectors, i + this->NbrProjectors - this->InitialNbrProjectors) = this->DiagonalizedMatrix(i, i);
	}
      return this->FullDiagonalizedMatrix;
    }
}

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::WriteState()
{
  if (this->Index > 1)
    system("mv lanczos.dat lanczos.bak");
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->AutomaticProjectorConstructionFlag);
  WriteLittleEndian(File, this->NbrProjectors);
  WriteLittleEndian(File, this->InitialNbrProjectors);
  if (this->ProjectorEigenvalues != 0)
    {
      for (int i = 0; i < (this->NbrProjectors - this->InitialNbrProjectors); ++i)
	WriteLittleEndian(File, this->ProjectorEigenvalues[i]);
    }
  int TmpDimension = this->TridiagonalizedMatrix.GetNbrRow();
  WriteLittleEndian(File, TmpDimension);
  --TmpDimension;
  for (int i = 0; i <= TmpDimension; ++i)    
    {    
      WriteLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
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

bool ComplexBasicLanczosAlgorithmWithGroundStateAndProjectorFastDisk::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->AutomaticProjectorConstructionFlag);
  ReadLittleEndian(File, this->NbrProjectors);
  ReadLittleEndian(File, this->InitialNbrProjectors);
  if (this->ProjectorEigenvalues != 0)
    {
      for (int i = 0; i < this->NbrProjectors; ++i)
	ReadLittleEndian(File, this->ProjectorEigenvalues[i]);
    }
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
  --TmpDimension;
  for (int i = 0; i <= TmpDimension; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  for (int i = 0; i < TmpDimension; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  File.close();
  char* TmpVectorName = new char [256];
  sprintf(TmpVectorName, "vector.%d", this->Index);
  this->V1.ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 1));
  this->V2.ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 2));
  this->V3.ReadVector(TmpVectorName);
  for (int i = 0; i < this->NbrProjectors; ++i)
    {
      sprintf(TmpVectorName, "vectorprojector.%d", i);
      this->ProjectorVectors[i].ReadVector(TmpVectorName);
    }
  delete[] TmpVectorName;
  return true;
}
