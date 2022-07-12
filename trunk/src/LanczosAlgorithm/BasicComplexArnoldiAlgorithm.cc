////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 17/11/2012                      //
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


#include "LanczosAlgorithm/BasicComplexArnoldiAlgorithm.h"
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

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
// leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

BasicComplexArnoldiAlgorithm::BasicComplexArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter, 
							   bool highEnergy, bool leftFlag, bool strongConvergence, bool sortRealFlag)  
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->ArnoldiVectors = new ComplexVector [this->MaximumNumberIteration];
  this->TemporaryCoefficients = new Complex [this->MaximumNumberIteration];
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = ComplexUpperHessenbergMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix();
      this->ReducedMatrix = ComplexUpperHessenbergMatrix();
   }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->StrongConvergenceFlag = strongConvergence;
  this->SortEigenvalueRealPartFlag = sortRealFlag;
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

BasicComplexArnoldiAlgorithm::BasicComplexArnoldiAlgorithm(const BasicComplexArnoldiAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->ArnoldiVectors = new ComplexVector [this->MaximumNumberIteration];
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

BasicComplexArnoldiAlgorithm::~BasicComplexArnoldiAlgorithm() 
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

void BasicComplexArnoldiAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectors[0] = ComplexVector (Dimension);
  this->ArnoldiVectors[1] = ComplexVector (Dimension);
  this->ArnoldiVectors[2] = ComplexVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    {
      this->ArnoldiVectors[0][i] = drand48() * Phase (2.0 * M_PI * drand48());
    }
  this->ArnoldiVectors[0] /= this->ArnoldiVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicComplexArnoldiAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectors[0] = vector;
  this->ArnoldiVectors[1] = ComplexVector (Dimension);
  this->ArnoldiVectors[2] = ComplexVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& BasicComplexArnoldiAlgorithm::GetGroundState()
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
  cout << "lapack is required for BasicComplexArnoldiAlgorithm" << endl;
#endif
  if (this->HighEnergyFlag == false)
    SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector, true);
  else
    SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector, true);

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

Vector* BasicComplexArnoldiAlgorithm::GetEigenstates(int nbrEigenstates)
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
  cout << "lapack is required for BasicComplexArnoldiAlgorithm" << endl;
#endif
  if (this->SortEigenvalueRealPartFlag == false)
    {
      if (this->HighEnergyFlag == false)
	SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector, true, 1e-10);
      else
	SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector, true, 1e-10);
    }
  else
    {
      if (this->HighEnergyFlag == false)
	SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector, false, 1e-10);
      else
	SortedDiagonalizedMatrix.SortMatrixDownOrder(TmpEigenvector, false, 1e-10);
    }

  Complex* TmpCoefficents = new Complex [SortedDiagonalizedMatrix.GetNbrColumn()];
  for (int i = 0; i < nbrEigenstates; ++i)
    {
      double TmpNorm = 1.0 / TmpEigenvector[i].Norm();
      for (int j = 0; j < SortedDiagonalizedMatrix.GetNbrColumn(); ++j)
	TmpCoefficents[j] = TmpEigenvector[i][j]*TmpNorm;
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

void BasicComplexArnoldiAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      if (nbrIter < 3)
	nbrIter = 3;
      Dimension = nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);

      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->ArnoldiVectors[0]), &(this->ArnoldiVectors[1]));
      Operation1.ApplyOperation(this->Architecture);
      this->ReducedMatrix.SetMatrixElement(0, 0, (this->ArnoldiVectors[0] * this->ArnoldiVectors[1]));
      Complex Tmp;
      this->ReducedMatrix.GetMatrixElement(0, 0, Tmp);
      this->ArnoldiVectors[1].AddLinearCombination(-Tmp, this->ArnoldiVectors[0]);
      this->ReducedMatrix.SetMatrixElement(1, 0, this->ArnoldiVectors[1].Norm()); 
      this->ReducedMatrix.GetMatrixElement(1, 0, Tmp);
      this->ArnoldiVectors[1] /=  Tmp; 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->ArnoldiVectors[1]), &(this->ArnoldiVectors[2]));
      Operation2.ApplyOperation(this->Architecture);
      this->ReducedMatrix.SetMatrixElement(0, 1, (this->ArnoldiVectors[0] * this->ArnoldiVectors[2]));
      this->ReducedMatrix.SetMatrixElement(1, 1, (this->ArnoldiVectors[1] * this->ArnoldiVectors[2]));
    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; ++i)
    {
      for (int k = 0; k < i; ++k)
	{
	  this->ReducedMatrix.GetMatrixElement(k, i - 1, this->TemporaryCoefficients[k]);
	  this->TemporaryCoefficients[k] *= -1.0;
	}
      AddComplexLinearCombinationOperation Operation2 (&(this->ArnoldiVectors[i]), this->ArnoldiVectors, i, this->TemporaryCoefficients);	  
      Operation2.ApplyOperation(this->Architecture);
      double VectorNorm = this->ArnoldiVectors[i].Norm();
      this->ReducedMatrix.SetMatrixElement(i, i - 1, VectorNorm);
      if (VectorNorm < 1e-5)
	{
	  cout << "subspace !!! " << i << endl;
	}
      this->ArnoldiVectors[i] /= VectorNorm;
      this->Index++;
      this->ArnoldiVectors[i + 1] = ComplexVector(this->Hamiltonian->GetHilbertSpaceDimension());
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->ArnoldiVectors[i]), &(this->ArnoldiVectors[i + 1]));
      Operation.ApplyOperation(this->Architecture);
      MultipleComplexScalarProductOperation Operation3 (&(this->ArnoldiVectors[i + 1]), this->ArnoldiVectors, i + 1, this->TemporaryCoefficients);
      Operation3.ApplyOperation(this->Architecture);
      for (int j = 0; j <= i; ++j)
	{
	  this->ReducedMatrix.SetMatrixElement(j, i, Conj(this->TemporaryCoefficients[j]));
	}
    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->ComplexPreviousWantedEigenvalues[i] = this->ComplexDiagonalizedMatrix[i];
      this->Diagonalize();
      if (this->SortEigenvalueRealPartFlag == false)
	{
	  if (this->HighEnergyFlag == false)
	    this->ComplexDiagonalizedMatrix.SortMatrixUpOrder(true, 1e-10);
	  else
	    this->ComplexDiagonalizedMatrix.SortMatrixDownOrder(true, 1e-10);
	}
      else
	{
	  if (this->HighEnergyFlag == false)
	    this->ComplexDiagonalizedMatrix.SortMatrixUpOrder();
	  else
	    this->ComplexDiagonalizedMatrix.SortMatrixDownOrder();
	}
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

bool BasicComplexArnoldiAlgorithm::TestConvergence ()
{
  if (this->ReducedMatrix.GetNbrRow() > this->NbrEigenvalue)
    {
      cout << this->Index << " : ";
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	cout << this->ComplexPreviousWantedEigenvalues[i] << " ";
      cout << Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1] - this->ComplexPreviousWantedEigenvalues[this->NbrEigenvalue - 1]) << endl;
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

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void BasicComplexArnoldiAlgorithm::GetEigenvalues (double*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new double [nbrEigenvalues];
  for (int i = 0; i < nbrEigenvalues; ++i)
    {
      eigenvalues[i] = this->ComplexDiagonalizedMatrix[i].Re;
    }
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void BasicComplexArnoldiAlgorithm::GetEigenvalues (Complex*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new Complex [nbrEigenvalues];
  for (int i = 0; i < nbrEigenvalues; ++i)
    {
      eigenvalues[i] = this->ComplexDiagonalizedMatrix[i];
    }
}

// diagonalize tridiagonalized matrix and find ground state energy
//

void BasicComplexArnoldiAlgorithm::Diagonalize () 
{
  int Dimension = this->ReducedMatrix.GetNbrRow();
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (this->TemporaryReducedMatrix.GetNbrColumn());
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag);
  this->ComplexDiagonalizedMatrix.Resize(this->TemporaryReducedMatrix.GetNbrColumn(), this->TemporaryReducedMatrix.GetNbrColumn());
  for (int i = 0; i < this->TemporaryReducedMatrix.GetNbrColumn(); ++i)
    this->ComplexDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "error, LAPACK is required for BasicComplexArnoldiAlgorithm" << endl;
#endif
  this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[0]);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (Norm(this->ComplexDiagonalizedMatrix[DiagPos]) < this->GroundStateEnergy)
      this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[DiagPos]);  
  return;
}

// restart the Arnoldi algorithm if needed
//

void BasicComplexArnoldiAlgorithm::RestartAlgorithm()
{
}
