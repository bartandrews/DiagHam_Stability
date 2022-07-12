////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of basic Lanczos algorithm                     //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 30/04/2001                      //
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


#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

FullReorthogonalizedLanczosAlgorithm::FullReorthogonalizedLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter,
									   bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->PreviousLastWantedEigenvalue = 0.0;
  this->StrongConvergenceFlag = strongConvergence;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

FullReorthogonalizedLanczosAlgorithm::FullReorthogonalizedLanczosAlgorithm(const FullReorthogonalizedLanczosAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
}

// destructor
//

FullReorthogonalizedLanczosAlgorithm::~FullReorthogonalizedLanczosAlgorithm() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
    }
  delete[] this->PreviousWantedEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void FullReorthogonalizedLanczosAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = RealVector (Dimension);
  this->LanczosVectors[1] = RealVector (Dimension);
  this->LanczosVectors[2] = RealVector (Dimension);
  int Shift = RAND_MAX / 2;
  double Scale = 1.0 / ((double) Shift);
  for (int i = 0; i < Dimension; i++)
    {
      this->LanczosVectors[0][i] = Scale * ((double) (rand() - Shift));
    }
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void FullReorthogonalizedLanczosAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = vector;
  this->LanczosVectors[1] = RealVector (Dimension);
  this->LanczosVectors[2] = RealVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& FullReorthogonalizedLanczosAlgorithm::GetGroundState()
{
  RealVector TmpComponents (this->DiagonalizedMatrix.GetNbrRow());
  this->TridiagonalizedMatrix.Eigenvector(this->GroundStateEnergy, TmpComponents);
/*  RealVector TmpComponents2(TmpComponents, true);
  TmpComponents2 *= this->TridiagonalizedMatrix;
  TmpComponents2 /= this->GroundStateEnergy;
  for (int i = 0; i < this->DiagonalizedMatrix.GetNbrRow(); i++)
    if (TmpComponents2[i] != TmpComponents[i])
      cout << i << " : " << TmpComponents2[i] << " " << TmpComponents[i] << endl;*/ 
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

Vector* FullReorthogonalizedLanczosAlgorithm::GetEigenstates(int nbrEigenstates)
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

void FullReorthogonalizedLanczosAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->LanczosVectors[0]), &(this->LanczosVectors[1]));
      Operation1.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->LanczosVectors[0] * 
							    this->LanczosVectors[1]);
      this->LanczosVectors[1].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index), 
						   this->LanczosVectors[0]);
      this->LanczosVectors[1] /= this->LanczosVectors[1].Norm(); 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[1]), &(this->LanczosVectors[2]));
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[0] * 
								       this->LanczosVectors[2]);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[1] * 
								      this->LanczosVectors[2]);
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      this->LanczosVectors[i].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index + 1), 
						   this->LanczosVectors[i - 1], 
						   -this->TridiagonalizedMatrix.
						   UpperDiagonalElement(this->Index), 
						   this->LanczosVectors[i - 2]);
      if (i > 2)
	{
	  double* TmpCoef = new double [i];
	  MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[i]), this->LanczosVectors, i - 2, TmpCoef);
	  Operation4.ApplyOperation(this->Architecture);
	  for (int j = 0; j < (i - 2); j++)
	    {
	      TmpCoef[j] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[i]), this->LanczosVectors, i - 2, TmpCoef);
	  Operation2.ApplyOperation(this->Architecture);
	  delete[] TmpCoef;
	}
      double VectorNorm = this->LanczosVectors[i].Norm();
      while (VectorNorm < 1e-5)
	{
	  cout << "subspace !!! " << i << endl;
	  double tmp = 0;
	  for (int j = 0; j < this->LanczosVectors[0].GetVectorDimension(); j++)
	    {
	      this->LanczosVectors[i][j] = (rand () - 16384) * 0.5;
	      tmp += this->LanczosVectors[i][j] * this->LanczosVectors[i][j];
	    }
	  tmp = sqrt(tmp);
	  this->LanczosVectors[i] /= tmp;
	  RealVector TmpVector(this->LanczosVectors[0].GetVectorDimension());
	  this->Hamiltonian->Multiply(this->LanczosVectors[i], TmpVector);
	  this->LanczosVectors[i] = TmpVector;

	  double* TmpCoef2 = new double [i];
	  for (int j = 0; j < i; j++)
	    {
	      TmpCoef2[j] = -(this->LanczosVectors[j] * this->LanczosVectors[i]);
	    }
	  AddRealLinearCombinationOperation Operation3 (&(this->LanczosVectors[i]), this->LanczosVectors, i, TmpCoef2);
	  Operation3.ApplyOperation(this->Architecture);
	  delete[] TmpCoef2;

/*	  for (int j = 0; j < i; j++)
	    {
	      this->LanczosVectors[i].AddLinearCombination(- (this->LanczosVectors[j] * this->LanczosVectors[i]), this->LanczosVectors[j]);
	    }*/
	  VectorNorm = this->LanczosVectors[i].Norm();
	}
      this->LanczosVectors[i] /= VectorNorm;
      this->Index++;
      this->LanczosVectors[i + 1] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[i]), &(this->LanczosVectors[i + 1]));
      Operation.ApplyOperation(this->Architecture);

      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[i - 1] * 
								       this->LanczosVectors[i + 1]);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[i] * 
								      this->LanczosVectors[i + 1]);

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

  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool FullReorthogonalizedLanczosAlgorithm::TestConvergence ()
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
		  if (fabs(this->DiagonalizedMatrix.DiagonalElement(i))>3*MACHINE_PRECISION)
		    return false;
		  else
		    if (fabs(this->PreviousWantedEigenvalues[i])>3*MACHINE_PRECISION)
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

