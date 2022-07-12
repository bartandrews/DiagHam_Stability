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


#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundState.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"

#include "GeneralTools/Endian.h"
#include "GeneralTools/ConfigurationParser.h"

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
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)

BasicLanczosAlgorithmWithGroundState::BasicLanczosAlgorithmWithGroundState(AbstractArchitecture* architecture, int maxIter, bool diskFlag, bool resumeDiskFlag) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->V1 = RealVector();
  this->V2 = RealVector();
  this->V3 = RealVector();
  this->InitialState = RealVector();
  this->GroundState = RealVector();
  this->GroundStateFlag = false;
  this->DiskFlag = diskFlag;
  this->ResumeDiskFlag = resumeDiskFlag;
  this->OrthogonalizationSetSize = 0;
  this->OrthogonalizationSet = 0;
  this->OrthogonalizationSetFileNames = 0;
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

// copy constructor
//
// algorithm = algorithm from which new one will be created

BasicLanczosAlgorithmWithGroundState::BasicLanczosAlgorithmWithGroundState(const BasicLanczosAlgorithmWithGroundState& algorithm) 
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
  this->NbrEigenvalue = 1;
  this->OrthogonalizationSetSize = 0;
  this->OrthogonalizationSet = 0;
  this->OrthogonalizationSetFileNames = 0;
}

// destructor
//

BasicLanczosAlgorithmWithGroundState::~BasicLanczosAlgorithmWithGroundState() 
{
  if (this->OrthogonalizationSet != 0)
    delete[] this->OrthogonalizationSet;
  if (this->OrthogonalizationSetFileNames != 0)
    {
      for (int i = 0; i < this->OrthogonalizationSetSize; ++i)
	delete[] this->OrthogonalizationSetFileNames[i];
      delete[] this->OrthogonalizationSetFileNames;
    }

}

// initialize Lanczos algorithm with a random vector
//

void BasicLanczosAlgorithmWithGroundState::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = RealVector (Dimension);
  this->V2 = RealVector (Dimension);
  this->V3 = RealVector (Dimension);
  if (this->ResumeDiskFlag == false)
    {
      int Shift = RAND_MAX / 2;
      double Scale = 1.0 / ((double) Shift);
      for (int i = 0; i < Dimension; i++)
	{
	  this->V1[i] = Scale * ((double) (rand() - Shift));
	}
      this->ExternalOrthonogalization(this->V1);
      this->V1 /= this->V1.Norm();
      if (this->DiskFlag == false)
	this->InitialState = RealVector (this->V1, true);
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

void BasicLanczosAlgorithmWithGroundState::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if (this->ResumeDiskFlag == false)
    {
      this->V1 = vector;
      this->V2 = RealVector (Dimension);
      this->V3 = RealVector (Dimension);
      if (this->OrthogonalizationSetSize > 0)
	{
	  this->ExternalOrthonogalization(this->V1);
	  this->V1 /= this->V1.Norm();
	}
      if (this->DiskFlag == false)
	this->InitialState = RealVector (vector, true);
      else
	this->V1.WriteVector("vector.0");
      this->Index = 0;
      this->GroundStateFlag = false;
      this->TridiagonalizedMatrix.Resize(0, 0);
    }
  else
    {
      this->V1 = RealVector (Dimension);
      this->V2 = RealVector (Dimension);
      this->V3 = RealVector (Dimension);
      this->ReadState();
    }
}

// force orthogonalization with respect to a set of vectors
//
// fileName = name of the file describing the set of vectors
// return value = true if no error occured

bool BasicLanczosAlgorithmWithGroundState::ForceOrthogonalization(char* fileName)
{
  ConfigurationParser OrthogonalizationSet;
  if (OrthogonalizationSet.Parse(fileName) == false)
    {
      this->OrthogonalizationSetSize = 0;
      OrthogonalizationSet.DumpErrors(cout) << endl;
      return false;
    }
  if (OrthogonalizationSet.GetAsStringArray("Vectors", ' ', this->OrthogonalizationSetFileNames, this->OrthogonalizationSetSize) == false)
    {
      cout << "Vectors are not defined or have a wrong value in " << fileName << endl;
      return false;
    }
  if (this->DiskFlag == false)
    {
      this->OrthogonalizationSet = new RealVector[this->OrthogonalizationSetSize];
      for (int i = 0; i < this->OrthogonalizationSetSize; ++i)
	this->OrthogonalizationSet[i].ReadVector(this->OrthogonalizationSetFileNames[i]);
    }
  return true;
}

// get the n first eigenstates (limited to the ground state fro this class, return NULL if nbrEigenstates > 1)
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* BasicLanczosAlgorithmWithGroundState::GetEigenstates(int nbrEigenstates)
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

Vector& BasicLanczosAlgorithmWithGroundState::GetGroundState()
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

      if (this->DiskFlag == false)
	{
	  double* TmpCoefficient = new double[2];
	  this->GroundState.Copy(this->InitialState, TmpComponents[0]);
	  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->InitialState, &this->V3);
	  Operation1.ApplyOperation(this->Architecture);
	  this->V3.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(0), this->InitialState);
	  this->ExternalOrthonogalization(this->V3);
	  this->V3 /= this->V3.Norm();
	  this->V2.Copy(this->InitialState);
	  this->GroundState.AddLinearCombination(TmpComponents[1], this->V3);
	  for (int i = 2; i < this->DiagonalizedMatrix.GetNbrRow(); ++i)
	    {
	      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V3, &this->V1);
	      Operation1.ApplyOperation(this->Architecture);
  	      RealVector* TmpVector = new RealVector[2];
	      TmpVector[0] = this->V2;
	      TmpVector[1] = this->V3;
	      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(i - 2);
	      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(i - 1);
	      AddRealLinearCombinationOperation Operation4 (&(this->V1),  TmpVector, 2, TmpCoefficient);
	      Operation4.ApplyOperation(this->Architecture);
	      delete[] TmpVector;
	      this->ExternalOrthonogalization(this->V1);
	      this->V1 /= this->V1.Norm();
	      this->GroundState.AddLinearCombination(TmpComponents[i], this->V1);
	      RealVector TmpV (this->V2);
	      this->V2 = this->V3;
	      this->V3 = this->V1;
	      this->V1 = TmpV;
 	      cout << i << "/" << this->DiagonalizedMatrix.GetNbrRow() << "           \r";
 	      cout.flush();
	    }
	}
      else
	{ 
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
      this->ExternalOrthonogalization(this->GroundState);
      this->GroundState /= this->GroundState.Norm();
      this->GroundStateFlag = true;
      delete[] TmpComponents;
    }
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicLanczosAlgorithmWithGroundState::RunLanczosAlgorithm (int nbrIter) 
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
      this->ExternalOrthonogalization(this->V2);
      this->V2 /= this->V2.Norm(); 
      if (this->DiskFlag == true)
	this->V2.WriteVector("vector.1");
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
  double* TmpCoefficient = new double[2];
  double TmpScalarProduct[2];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      RealVector* TmpVector = new RealVector[2];
      if (this->ResumeDiskFlag == false)
	{
	  TmpVector[0] = this->V1;
	  TmpVector[1] = this->V2;
	  TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
	  TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
	  AddRealLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
	  Operation4.ApplyOperation(this->Architecture);
	  delete[] TmpVector;
	  this->ExternalOrthonogalization(this->V3);
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
	  RealVector TmpV (this->V2);
	  this->V2 = this->V3;
	  this->V3 = TmpV;	  
	  this->V1 = RealVector();
	}
      else
	{
	  RealVector TmpV (this->V1);
	  this->V1 = this->V2;
	  this->V2 = this->V3;
	  this->V3 = TmpV;
	}
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V2, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
      if (this->DiskFlag == true)
	{
	  char* TmpVectorName = new char [256];
	  sprintf(TmpVectorName, "vector.%d", (i - 1));
	  this->V1.ReadVector(TmpVectorName);
	  delete[] TmpVectorName;
	}      
      RealVector* TmpVectorScalarProduct[2];
      TmpVectorScalarProduct[0] = &(this->V1);
      TmpVectorScalarProduct[1] = &(this->V2);
      MultipleRealScalarProductOperation Operation2 (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = TmpScalarProduct[0];
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = TmpScalarProduct[1];
    }
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
}
  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool BasicLanczosAlgorithmWithGroundState::TestConvergence ()
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

bool BasicLanczosAlgorithmWithGroundState::WriteState()
{
  ofstream File;
  File.open("lanczos.dat", ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
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

bool BasicLanczosAlgorithmWithGroundState::ReadState()
{
  ifstream File;
  File.open("lanczos.dat", ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
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
  delete[] TmpVectorName;
  return true;
}


// orthogonalize a vector with respect to a set of external vectors
//
// inputVector = reference on the vector whose component on the external set has to be removed 

void  BasicLanczosAlgorithmWithGroundState::ExternalOrthonogalization(RealVector& inputVector)
{
  if (this->OrthogonalizationSetSize > 0)
    {
      if (this->DiskFlag == false)
	{
	  double* TmpCoefficient = new double[this->OrthogonalizationSetSize];
	  for (int i = 0; i < this->OrthogonalizationSetSize; ++i)
	    TmpCoefficient[i] = - (inputVector * this->OrthogonalizationSet[i]);
	  AddRealLinearCombinationOperation Operation (&inputVector, this->OrthogonalizationSet, 
						       this->OrthogonalizationSetSize, TmpCoefficient);
	  Operation.ApplyOperation(this->Architecture);	      
	  delete[] TmpCoefficient;
	}
      else
	{
	  double TmpCoefficient = 0.0;
	  RealVector TmpVector(this->Hamiltonian->GetHilbertSpaceDimension());
	  for (int i = 0; i < this->OrthogonalizationSetSize; ++i)
	    {
	      TmpVector.ReadVector(this->OrthogonalizationSetFileNames[i]);	      
	      TmpCoefficient = - (inputVector * TmpVector);
	      AddRealLinearCombinationOperation Operation (&inputVector, &(TmpVector), 
							   1, &TmpCoefficient);
	      Operation.ApplyOperation(this->Architecture);	      
	    }
	}
    }
}
