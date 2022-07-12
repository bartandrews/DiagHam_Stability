////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of full reorthogonalized Lanczos algorithm             //
//                (with full re-orthogonalization at each step)               //
//                 and storing each iteration information on disk             //
//                                                                            //
//                        last modification : 18/03/2003                      //
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


#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


using std::cout;
using std::endl;
using std::ios;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

FullReorthogonalizedLanczosAlgorithmWithDiskStorage::FullReorthogonalizedLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, 
													 int nbrEigenvalue, int maxNbrVectors, int maxIter,
													 bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->MaxNbrVectors = maxNbrVectors;
  this->LanczosVectors = new RealVector [this->MaxNbrVectors];
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
  this->StrongConvergenceFlag = strongConvergence;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

FullReorthogonalizedLanczosAlgorithmWithDiskStorage::FullReorthogonalizedLanczosAlgorithmWithDiskStorage(const FullReorthogonalizedLanczosAlgorithmWithDiskStorage& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->MaxNbrVectors = algorithm.MaxNbrVectors;
  this->LanczosVectors = new RealVector [this->MaxNbrVectors];
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

FullReorthogonalizedLanczosAlgorithmWithDiskStorage::~FullReorthogonalizedLanczosAlgorithmWithDiskStorage() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
    }
}

// initialize Lanczos algorithm with a random vector
//

void FullReorthogonalizedLanczosAlgorithmWithDiskStorage::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->LanczosVectors[i] = RealVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    this->LanczosVectors[0][i] = (rand() - 32767) * 0.5;
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void FullReorthogonalizedLanczosAlgorithmWithDiskStorage::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->LanczosVectors[i] = RealVector (Dimension);
  this->LanczosVectors[0] = vector;
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// resume Lanczos algorithm from disk datas in current directory
//

void FullReorthogonalizedLanczosAlgorithmWithDiskStorage::ResumeLanczosAlgorithm()
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->LanczosVectors[i] = RealVector (Dimension);
  this->ReadState();
  char* TmpVectorName = new char [256];
  
  sprintf(TmpVectorName, "vector.%d", this->Index);
  this->LanczosVectors[0].ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 1));
  this->LanczosVectors[1].ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index + 2));
  this->LanczosVectors[2].ReadVector(TmpVectorName);
  sprintf(TmpVectorName, "vector.%d", (this->Index - 1));
  this->LanczosVectors[3].ReadVector(TmpVectorName);
}
  
// get last produced vector
//
// return value = reference on lest produced vector

Vector& FullReorthogonalizedLanczosAlgorithmWithDiskStorage::GetGroundState()
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

Vector* FullReorthogonalizedLanczosAlgorithmWithDiskStorage::GetEigenstates(int nbrEigenstates)
{
  this->Index = 0;
  RealVector* Eigenstates = new RealVector [nbrEigenstates];
  RealMatrix TmpEigenvector (this->TridiagonalizedMatrix.GetNbrRow(), this->TridiagonalizedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->TridiagonalizedMatrix.GetNbrRow());
  SortedDiagonalizedMatrix.Copy(this->TridiagonalizedMatrix);
  SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  double* TmpCoefficents = new double [this->TridiagonalizedMatrix.GetNbrRow()];
  char* TmpVectorName = new char [256];
  this->LanczosVectors[0].ReadVector("vector.0");
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      for (int j = 0; j < this->TridiagonalizedMatrix.GetNbrRow(); ++j)
	TmpCoefficents[j] = TmpEigenvector(j, i);
      Eigenstates[i] = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
      Eigenstates[i].Copy(this->LanczosVectors[0], TmpEigenvector(0, i));
      int ReducedMaxNbrVector =  this->MaxNbrVectors - 1;
      int MaxPos = (this->TridiagonalizedMatrix.GetNbrRow() - 1) / ReducedMaxNbrVector;
      int k = 0;
      for (; k < MaxPos; ++k)
	{
	  for (int j = 0; j < ReducedMaxNbrVector; ++j)
	    {
	      sprintf(TmpVectorName, "vector.%d", (1 + j + (k * ReducedMaxNbrVector)));
	      this->LanczosVectors[1 + j].ReadVector(TmpVectorName);
	    }
	  AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), ReducedMaxNbrVector, &(TmpCoefficents[1 + (k * ReducedMaxNbrVector)]));
	  Operation.ApplyOperation(this->Architecture);
	}
      MaxPos = (this->TridiagonalizedMatrix.GetNbrRow() - 1) - (MaxPos * ReducedMaxNbrVector);
      for (int j = 0; j < MaxPos; ++j)
	{
	  sprintf(TmpVectorName, "vector.%d", (1 + j + (k * ReducedMaxNbrVector)));
	  this->LanczosVectors[1 + j].ReadVector(TmpVectorName);
	}
      AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), MaxPos, &(TmpCoefficents[1 + (k * ReducedMaxNbrVector)]));
      Operation.ApplyOperation(this->Architecture);
      Operation.ApplyOperation(this->Architecture);
      Eigenstates[i] /= Eigenstates[i].Norm();
    }
  delete[] TmpVectorName;
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void FullReorthogonalizedLanczosAlgorithmWithDiskStorage::RunLanczosAlgorithm (int nbrIter) 
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
      this->LanczosVectors[0].WriteVector("vector.0");
      this->LanczosVectors[1].WriteVector("vector.1");
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  char* TmpVectorName = new char [256];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      this->LanczosVectors[2].AddLinearCombination(-this->TridiagonalizedMatrix.
						   DiagonalElement(this->Index + 1), 
						   this->LanczosVectors[1], 
						   -this->TridiagonalizedMatrix.
						   UpperDiagonalElement(this->Index), 
						   this->LanczosVectors[0]);
      if (i > 2)
	{
	  double* TmpCoef = new double [i];
	  int ReducedMaxNbrVector = this->MaxNbrVectors - 4;
	  int MaxPos = (i - 3) / ReducedMaxNbrVector;
	  int k = 0;
	  for (; k < MaxPos; ++k)
	    {
	      for (int j = 0; j < ReducedMaxNbrVector; ++j)
		{
		  sprintf(TmpVectorName, "vector.%d", (j + (k * ReducedMaxNbrVector)));
		  this->LanczosVectors[4 + j].ReadVector(TmpVectorName);
		}
	      MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[2]), &(this->LanczosVectors[4]), ReducedMaxNbrVector, TmpCoef);
	      Operation4.ApplyOperation(this->Architecture);
	      for (int j = 0; j < ReducedMaxNbrVector; j++)
		{
		  TmpCoef[j] *= -1.0;
		}
	      AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[2]), &(this->LanczosVectors[4]), ReducedMaxNbrVector, TmpCoef);
	      Operation2.ApplyOperation(this->Architecture);
	    }
	  MaxPos = (i - 3) - (MaxPos * ReducedMaxNbrVector);
	  for (int j = 0; j < MaxPos; ++j)
	    {
	      sprintf(TmpVectorName, "vector.%d", (j + (k * ReducedMaxNbrVector)));
	      this->LanczosVectors[4 + j].ReadVector(TmpVectorName);
	    }
	  MultipleRealScalarProductOperation Operation4 (&(this->LanczosVectors[2]), &(this->LanczosVectors[3]), MaxPos + 1, TmpCoef);
	  Operation4.ApplyOperation(this->Architecture);
	  for (int j = 0; j <= MaxPos; j++)
	    {
	      TmpCoef[j] *= -1.0;
	    }
	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[2]), &(this->LanczosVectors[3]), MaxPos + 1, TmpCoef);
	  Operation2.ApplyOperation(this->Architecture);	  
	  delete[] TmpCoef;
	}

      double VectorNorm = this->LanczosVectors[2].Norm();
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
	  this->Hamiltonian->Multiply(this->LanczosVectors[2], TmpVector);
	  this->LanczosVectors[i] = TmpVector;

	  double* TmpCoef2 = new double [i];
	  for (int j = 0; j < i; j++)
	    {
	      TmpCoef2[j] = -(this->LanczosVectors[j] * this->LanczosVectors[i]);
	    }
	  AddRealLinearCombinationOperation Operation3 (&(this->LanczosVectors[2]), this->LanczosVectors, i, TmpCoef2);
	  Operation3.ApplyOperation(this->Architecture);
	  delete[] TmpCoef2;

	  VectorNorm = this->LanczosVectors[i].Norm();
	}
      this->LanczosVectors[2] /= VectorNorm;

      RealVector TmpV (this->LanczosVectors[0]);
      this->LanczosVectors[0] = this->LanczosVectors[1];
      this->LanczosVectors[1] = this->LanczosVectors[2];
      this->LanczosVectors[2] = this->LanczosVectors[3];
      this->LanczosVectors[3] = TmpV;
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[1]), &(this->LanczosVectors[2]));
      Operation.ApplyOperation(this->Architecture);

      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->LanczosVectors[0] * 
								       this->LanczosVectors[2]);
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->LanczosVectors[1] * 
								      this->LanczosVectors[2]);
      sprintf(TmpVectorName, "vector.%d", i);
      this->LanczosVectors[1].WriteVector(TmpVectorName);
      sprintf(TmpVectorName, "vector.%d", i + 1);
      this->LanczosVectors[2].WriteVector(TmpVectorName);
      this->WriteState();
    }
  delete[] TmpVectorName;
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
  this->WriteState();
}

  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool FullReorthogonalizedLanczosAlgorithmWithDiskStorage::TestConvergence ()
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

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool FullReorthogonalizedLanczosAlgorithmWithDiskStorage::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->EigenvaluePrecision);
  WriteLittleEndian(File, this->NbrEigenvalue);
  int TmpDimension = this->TridiagonalizedMatrix.GetNbrRow();
  WriteLittleEndian(File, TmpDimension);
  for (int i = 0; i <= (this->Index + 1); ++i)    
    {    
      WriteLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  for (int i = 0; i <= this->Index; ++i)
    {
      WriteLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
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

bool FullReorthogonalizedLanczosAlgorithmWithDiskStorage::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->EigenvaluePrecision);
  ReadLittleEndian(File, this->NbrEigenvalue);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
  for (int i = 0; i <= (this->Index + 1); ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  for (int i = 0; i <= this->Index; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  if (this->PreviousWantedEigenvalues != 0)
    delete[] this->PreviousWantedEigenvalues;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    {
      ReadLittleEndian(File, this->PreviousWantedEigenvalues[i]);
      this->PreviousWantedEigenvalues[i] *= 2.0;
    }
  File.close();  
  this->Diagonalize();
  this->DiagonalizedMatrix.SortMatrixUpOrder();
  return true;
}


