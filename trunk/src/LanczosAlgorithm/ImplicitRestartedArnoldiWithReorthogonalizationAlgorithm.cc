////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of implicit restarted Arnoldi algorithm             //
//                          (with re-orthogonalization)                       //
//                                                                            //
//                        last modification : 06/01/2003                      //
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


#include "LanczosAlgorithm/ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Matrix/RealMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvectors = number of needed eigenvectors (and eigenvalues)
// nbrIteration = number of Lanczos iteration per Arnoldi iteration (equal 2 * nbrEigenvectors)

ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm(AbstractArchitecture* architecture,
														   int nbrEigenvectors, int nbrIteration) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->NbrEigenvectors = nbrEigenvectors;
  if (nbrIteration != 0)
    this->NbrIteration = nbrIteration;
  else
    this->NbrIteration = 2 * nbrEigenvectors;
  if (this->NbrIteration > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->NbrIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->NbrIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->NbrUnwantedEigenvalues = this->NbrIteration - this->NbrEigenvectors;
  this->UnwantedEigenvalues =  new double [this->NbrUnwantedEigenvalues];
  this->Architecture = architecture;
  this->ConvergedValues = new double [this->NbrEigenvectors];
  this->NbrConvergedValue = 0;
  this->LastErrorRitzValue = 1.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm(const 
														   ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm& 
														   algorithm) 
{
  this->Index = algorithm.Index;
  this->NbrEigenvectors = algorithm.NbrEigenvectors;
  this->NbrIteration = algorithm.NbrIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->LanczosVectors = algorithm.LanczosVectors;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->NbrUnwantedEigenvalues = algorithm.NbrUnwantedEigenvalues;
  this->UnwantedEigenvalues =  new double [this->NbrUnwantedEigenvalues];
  this->ConvergedValues = new double [this->NbrEigenvectors];
  this->NbrConvergedValue = 0;
  this->LastErrorRitzValue = 0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// destructor
//

ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::~ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm() 
{
  delete[] this->UnwantedEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors = RealMatrix (Dimension, this->NbrIteration + 1, true);
  this->ConvergedEigenvectors  = RealMatrix (Dimension, this->NbrEigenvectors, true);
  for (int i = 0; i < Dimension; i++)
    this->LanczosVectors[0][i] = (rand() - 32767) * 0.5;
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  this->Index = 0;
//  this->StandardArnoldiStep();
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors = RealMatrix (Dimension, this->NbrIteration + 1, true);
  this->ConvergedEigenvectors  = RealMatrix (Dimension, this->NbrEigenvectors, true);
  this->LanczosVectors[0] = vector;
  this->Index = 0;
//  this->StandardArnoldiStep();
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::GetGroundState()
{
  RealVector TmpComponents (this->DiagonalizedMatrix.GetNbrRow());
  this->TridiagonalizedMatrix.Eigenvector(this->GroundStateEnergy, TmpComponents);
  this->GroundState.Copy(this->LanczosVectors[0], TmpComponents[0]);
  for (int i = 1; i < this->DiagonalizedMatrix.GetNbrRow(); i++)
    this->GroundState.AddLinearCombination (TmpComponents[i], this->LanczosVectors[i]);
  this->GroundState /= this->GroundState.Norm();
  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::GetEigenstates(int nbrEigenstates)
{
  RealVector* Eigenstates = new RealVector [nbrEigenstates];
  RealMatrix TmpEigenvector (this->TridiagonalizedMatrix.GetNbrRow(), this->TridiagonalizedMatrix.GetNbrRow(), true);  
  for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->TridiagonalizedMatrix.GetNbrRow());
  SortedDiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
  SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  double* TmpCoefficents = new double [this->TridiagonalizedMatrix.GetNbrRow()];
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      for (int j = 0; j < this->TridiagonalizedMatrix.GetNbrRow(); ++j)
	TmpCoefficents[j] = TmpEigenvector(j, i);
      Eigenstates[i] = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
      Eigenstates[i].Copy(this->LanczosVectors[0], TmpEigenvector(0, i));
      AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), this->TridiagonalizedMatrix.GetNbrRow() - 1, &(TmpCoefficents[1]));
      Operation.ApplyOperation(this->Architecture);
      Eigenstates[i] /= Eigenstates[i].Norm();
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
/*  for (int i = 0; i < 20; ++i)
    cout << this->LanczosVectors[0][i] << " ";
  cout << endl;*/
  cout << "new iteration.........." << endl;
  this->StandardArnoldiStep();

/*  cout << "norm = " << endl;
  for (int i = 0; i <= this->NbrIteration; ++i)
    cout << i << " = " << this->LanczosVectors[i].Norm() << endl;
  cout << "last norm = " << this->LastVectorNorm << endl;
  cout << "ortho = " << endl;
  for (int i = 0; i < this->NbrIteration; ++i)
    for (int j = i + 1; j <= this->NbrIteration; ++j)
      cout << i << " " << j << " = " << (this->LanczosVectors[i] * this->LanczosVectors[j]) << endl;*/

  // exact shift method
  this->DiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
  this->DiagonalizedMatrix.Diagonalize(100);
  this->DiagonalizedMatrix.SortMatrixUpOrder();
/*  for (int i = 0; i < this->NbrIteration; ++i)
    cout << this->DiagonalizedMatrix.DiagonalElement(i) << endl;*/
  for (int i = 0; i < this->NbrUnwantedEigenvalues; ++i)
    this->UnwantedEigenvalues[i] = this->DiagonalizedMatrix.DiagonalElement(this->TridiagonalizedMatrix.GetNbrRow() - 1 - i);

//  cout << this->TridiagonalizedMatrix << endl<< endl;

  // apply shifts  
  RealMatrix Q (this->NbrIteration, this->NbrIteration, true);
  for (int i = 0; i < this->NbrIteration; ++i)
    Q(i, i) = 1.0;
/*  for (int i = 0; i < (this->NbrIteration - 1); ++i)
    cout << i << " " << this->TridiagonalizedMatrix.DiagonalElement(i) << " " << this->TridiagonalizedMatrix(i, i + 1) << endl;
  cout << (this->NbrIteration - 1) << " " << this->TridiagonalizedMatrix.DiagonalElement(this->NbrIteration - 1) <<endl;*/
  RealTriDiagonalSymmetricMatrix TmpM (this->TridiagonalizedMatrix.PolynomialFilterWithExactShift(Q, this->UnwantedEigenvalues, this->NbrUnwantedEigenvalues));
//  cout << TmpM << endl;
/*  for (int i = 0; i < (this->NbrIteration - 1); ++i)
    cout << i << " " << TmpM.DiagonalElement(i) << " " << TmpM(i, i + 1) << endl;
  cout << (this->NbrIteration - 1) << " " << TmpM.DiagonalElement(this->NbrIteration - 1) <<endl;*/
/*  RealMatrix Q2 ((Matrix&) Q);
  Q2.Transpose();
  cout << (Q * Q2) << endl;*/
  this->TridiagonalizedMatrix = TmpM;
//  cout << Q << endl;
  // prepare restarting
  Q.Resize(this->NbrIteration, this->NbrEigenvectors + 1);
  this->LanczosVectors.Resize(this->Hamiltonian->GetHilbertSpaceDimension(), this->NbrIteration);
  MatrixMatrixMultiplyOperation Operation(&this->LanczosVectors, &Q);
  Operation.ApplyOperation(this->Architecture);
//  this->LanczosVectors.Multiply(Q);
  this->LanczosVectors.Resize(this->Hamiltonian->GetHilbertSpaceDimension(), this->NbrIteration + 1); 
  // construct last lanczos vector
  RealVector TmpLastVector = this->LanczosVectors[this->NbrIteration];
  this->LanczosVectors[this->NbrIteration] = this->LanczosVectors[this->NbrEigenvectors];
  this->LanczosVectors[this->NbrEigenvectors] = TmpLastVector;
/*  cout << "norm = " << this->LanczosVectors[this->NbrEigenvectors].Norm() << " " << this->LastVectorNorm << endl;
  cout << "Q factor = " << Q(this->NbrIteration - 1, this->NbrEigenvectors - 1) << endl;
  cout << "Q factor = " << Q(this->NbrIteration - 1, this->NbrEigenvectors - 2) << endl;
  cout << "Q factor = " << Q(this->NbrIteration - 1, this->NbrEigenvectors + 1) << endl;*/
  this->LastVectorNorm *= fabs(Q(this->TridiagonalizedMatrix.GetNbrRow() - 1, this->NbrEigenvectors - 1));
  this->LanczosVectors[this->NbrEigenvectors] *= Q(this->NbrIteration - 1, this->NbrEigenvectors - 1);
  double Factor = this->TridiagonalizedMatrix(this->NbrEigenvectors, this->NbrEigenvectors - 1);
  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
    this->LanczosVectors[this->NbrEigenvectors][i] += Factor * this->LanczosVectors[this->NbrIteration][i];
  this->LastVectorNorm  = this->LanczosVectors[this->NbrEigenvectors].Norm();
/*  for (int i = 0; i < this->NbrEigenvectors; ++i)
    {
      cout << (this->LanczosVectors[i] * this->LanczosVectors[this->NbrEigenvectors]) / this->LastVectorNorm << " ";
    }
  cout << endl;*/
  
  this->TridiagonalizedMatrix.Resize(this->NbrEigenvectors, this->NbrEigenvectors);
  RealMatrix Q2 (this->TridiagonalizedMatrix.GetNbrRow(), this->TridiagonalizedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
    Q2(i, i) = 1.0;
  this->DiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
  this->DiagonalizedMatrix.Diagonalize(Q2, 50);
  this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(0);
  // find error Ritz for each Ritz Value
  for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
    {
      cout << this->DiagonalizedMatrix.DiagonalElement(i) << " -> " << (fabs(Q2(this->TridiagonalizedMatrix.GetNbrRow() - 1, i)) * this->LastVectorNorm) << endl;
    }
  this->LastErrorRitzValue = (fabs(Q2(this->TridiagonalizedMatrix.GetNbrRow() - 1, this->TridiagonalizedMatrix.GetNbrRow() - 1)) * this->LastVectorNorm);
  bool TmpLockFlag = false;
  int TmpLockPos = 0;
  while ((TmpLockPos < this->TridiagonalizedMatrix.GetNbrRow()) && (TmpLockFlag == false))
    if ((fabs(Q2(this->TridiagonalizedMatrix.GetNbrRow() - 1,TmpLockPos)) * this->LastVectorNorm) < this->EigenvaluePrecision)
      TmpLockFlag = true;
    else
      ++TmpLockPos;
  if (TmpLockFlag == true)
    {
      this->ConvergedValues[this->NbrConvergedValue] = this->DiagonalizedMatrix.DiagonalElement(TmpLockPos);
      this->LanczosVectors.Resize(this->LanczosVectors.GetNbrRow(),  this->TridiagonalizedMatrix.GetNbrRow());
      MatrixMatrixMultiplyOperation LockFlagOperation(&this->LanczosVectors, &Q2);
      LockFlagOperation.ApplyOperation(this->Architecture);
      this->LanczosVectors.Resize(this->LanczosVectors.GetNbrRow(), this->NbrIteration + 1);
/*      for (int i = 0; i < TmpVector.GetVectorDimension(); ++i)
	{
	  TmpVector[i] = 0.0;
	  for (int j = 0; j < this->TridiagonalizedMatrix.GetNbrRow(); ++j)
	    TmpVector[i] += this->LanczosVectors[j][i] * Q2(j, TmpLockPos);
	}
      TmpVector /= TmpVector.Norm();  */
      RealVector TmpVector (this->ConvergedEigenvectors[this->NbrConvergedValue]);
      this->ConvergedEigenvectors[this->NbrConvergedValue] = this->LanczosVectors[TmpLockPos];
      this->LanczosVectors[TmpLockPos] = TmpVector;
      if (TmpLockPos == 0)
	{
	  TmpVector = this->LanczosVectors[1];
	  this->LanczosVectors[1] = this->LanczosVectors[TmpLockPos];
	  this->LanczosVectors[TmpLockPos] = this->LanczosVectors[1];
	}
      ++this->NbrConvergedValue;
//      for (int i = 0; i < TmpVector.GetVectorDimension(); ++i)
//	this->LanczosVectors[0][i] = (drand48() - 0.5) * 2.0;
      double* TmpCoef = new double [this->NbrConvergedValue];
      MultipleRealScalarProductOperation Operation4 (&this->LanczosVectors[0], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
      Operation4.ApplyOperation(this->Architecture);
      for (int j = 0; j < this->NbrConvergedValue; j++)
	{
	  TmpCoef[j] *= -1.0;
	}
      AddRealLinearCombinationOperation Operation2 (&this->LanczosVectors[0], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
      Operation2.ApplyOperation(this->Architecture);
      this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
      delete[] TmpCoef;
      this->Index = 0;
      cout << "value " << TmpLockPos << " has been locked" << endl;
    }
  else
    {
      this->Index = this->NbrEigenvectors;
    }
  this->TridiagonalizedMatrix.Resize(this->NbrIteration, this->NbrIteration);

  
/*  if (fabs(this->TridiagonalizedMatrix(0, 1)) < 1e-14)
    {
      cout << "found eigenvalue " << this->NbrConvergedValue << ": " << this->TridiagonalizedMatrix(0, 0) << endl;
      this->ConvergedValues[this->NbrConvergedValue] = this->TridiagonalizedMatrix(0, 0);
      RealVector TmpVector = this->ConvergedEigenvectors[this->NbrConvergedValue];
      this->ConvergedEigenvectors[this->NbrConvergedValue] = this->LanczosVectors[0];
      this->LanczosVectors[0] = TmpVector;
      ++this->NbrConvergedValue;
      for (int i = 0; i < TmpVector.GetVectorDimension(); ++i)
	TmpVector[i] = (drand48() - 0.5) * 2.0;
      double* TmpCoef = new double [this->NbrConvergedValue];
      MultipleRealScalarProductOperation Operation4 (&TmpVector, this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
      Operation4.ApplyOperation(this->Architecture);
      for (int j = 0; j < this->NbrConvergedValue; j++)
	{
	  TmpCoef[j] *= -1.0;
	}
      AddRealLinearCombinationOperation Operation2 (&TmpVector, this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
      Operation2.ApplyOperation(this->Architecture);
      TmpVector /= TmpVector.Norm();
      delete[] TmpCoef;
      this->Index = 0;
      --this->NbrEigenvectors;
      ++this->NbrUnwantedEigenvalues;
    }
  else*/
    {
      /*  cout << "norm = " << this->LanczosVectors[this->NbrEigenvectors].Norm() << " " << this->LastVectorNorm << endl;
	  cout << "norm = " << this->LanczosVectors[0].Norm() << " ortho = " <<( this->LanczosVectors[0] * this->LanczosVectors[1]) << endl;
	  cout << "norm = " << this->LanczosVectors[1].Norm()<< " ortho = " << (this->LanczosVectors[2] * this->LanczosVectors[1]) << endl;
	  cout << "norm = " << this->LanczosVectors[2].Norm()<< " ortho = " << (this->LanczosVectors[0] * this->LanczosVectors[2]) << endl;*/
/*      cout << "norm = " << endl;
      for (int i = 0; i <= this->NbrEigenvectors; ++i)
	cout << i << " = " << this->LanczosVectors[i].Norm() << endl;
      cout << "last norm = " << this->LastVectorNorm << endl;
      cout << "ortho = " << endl;
      for (int i = 0; i < this->NbrEigenvectors; ++i)
	for (int j = i + 1; j <= this->NbrEigenvectors; ++j)
	  cout << i << " " << j << " = " << (this->LanczosVectors[i] * this->LanczosVectors[j]) << endl;*/

//      this->Index = this->NbrEigenvectors;
//      this->Index = 0;
    }
  if (this->NbrConvergedValue > 0)
    {
      cout << "scalar = " << (this->ConvergedEigenvectors[0] * this->LanczosVectors[0]) << endl;
    }
}
  
// standard Arnoldi step with full-reothogonalization
//

void ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::StandardArnoldiStep()
{
  if (this->Index == 0)
    {
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->LanczosVectors[0]), &(this->LanczosVectors[1]));
      Operation1.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(0) = (this->LanczosVectors[0] * this->LanczosVectors[1]);
      this->LanczosVectors[1].AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(0), 
						   this->LanczosVectors[0]);

      if (this->NbrConvergedValue > 0)
	{
	  double* TmpCoef = new double [this->NbrConvergedValue];
	  MultipleRealScalarProductOperation Operation6 (&this->LanczosVectors[1], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
	  Operation6.ApplyOperation(this->Architecture);
	  for (int j = 0; j < this->NbrConvergedValue; j++)
	    {
	      TmpCoef[j] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation7 (&this->LanczosVectors[1], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
	  Operation7.ApplyOperation(this->Architecture);
	  delete[] TmpCoef;
	}

      this->LastVectorNorm = this->LanczosVectors[1].Norm();
      this->Index = 1;
    }
  if (this->Index == 1)
    {
      this->TridiagonalizedMatrix.UpperDiagonalElement(0) = this->LastVectorNorm;
      this->LanczosVectors[1] /= this->LastVectorNorm; 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[1]), &(this->LanczosVectors[2]));
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(1) = (this->LanczosVectors[1] * this->LanczosVectors[2]);
      this->LanczosVectors[2].AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(1), 
						   this->LanczosVectors[1], -this->TridiagonalizedMatrix.UpperDiagonalElement(0), 
						   this->LanczosVectors[0]);

      if (this->NbrConvergedValue > 0)
	{
	  double* TmpCoef = new double [this->NbrConvergedValue];
	  MultipleRealScalarProductOperation Operation5 (&this->LanczosVectors[2], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
	  Operation5.ApplyOperation(this->Architecture);
	  for (int j = 0; j < this->NbrConvergedValue; j++)
	    {
	      TmpCoef[j] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation3 (&this->LanczosVectors[2], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef);
	  Operation3.ApplyOperation(this->Architecture);
	  delete[] TmpCoef;
	}
      
      this->LastVectorNorm = this->LanczosVectors[2].Norm();
      this->Index = 2;

   }
  double* TmpCoef = new double [this->NbrIteration];
  for (int i = this->Index; i < this->NbrIteration; ++i)
    {
      this->TridiagonalizedMatrix.UpperDiagonalElement(i - 1) = this->LastVectorNorm;
      this->LanczosVectors[i] /= this->LastVectorNorm;
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[i]), &(this->LanczosVectors[i + 1]));
      Operation.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(i) = (this->LanczosVectors[i] * this->LanczosVectors[i + 1]);
      this->LanczosVectors[i + 1].AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(i), 
						       this->LanczosVectors[i], -this->TridiagonalizedMatrix.UpperDiagonalElement(i - 1), 
						       this->LanczosVectors[i - 1]);

      if (this->NbrConvergedValue > 0)
	{
	  double* TmpCoef2 = new double [this->NbrConvergedValue];
	  MultipleRealScalarProductOperation Operation6 (&this->LanczosVectors[i + 1], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef2);
	  Operation6.ApplyOperation(this->Architecture);
	  for (int j = 0; j < this->NbrConvergedValue; j++)
	    {
	      TmpCoef2[j] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation7 (&this->LanczosVectors[i + 1], this->ConvergedEigenvectors, this->NbrConvergedValue, TmpCoef2);
	  Operation7.ApplyOperation(this->Architecture);
	  delete[] TmpCoef2;
	}
      if (i > 2)
	{
//	  MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[i]), this->LanczosVectors, i - 1, TmpCoef);
	  MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[i + 1]), this->LanczosVectors, i + 1, TmpCoef);
	  Operation4.ApplyOperation(this->Architecture);
//	  for (int j = 0; j < (i - 1); j++)
	  for (int j = 0; j <= i; j++)
	    {
	      TmpCoef[j] *= -1.0;
	    }
//	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[i]), this->LanczosVectors, i - 1, TmpCoef);
	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[i + 1]), this->LanczosVectors, i + 1, TmpCoef);
	  Operation2.ApplyOperation(this->Architecture);
	}
      this->LastVectorNorm = this->LanczosVectors[i + 1].Norm();
      ++this->Index;
    } 
  delete[] TmpCoef;
}

// test if convergence has been reached
//
// return value = true if convergence has been reached

bool ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm::TestConvergence ()
{
  cout << this->LastErrorRitzValue << " " << this->EigenvaluePrecision << endl;
/*  if (this->LastErrorRitzValue  < this->EigenvaluePrecision)
    return true;
  else*/
    return false;
}
