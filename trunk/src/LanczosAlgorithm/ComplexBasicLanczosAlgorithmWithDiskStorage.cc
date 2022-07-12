////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of basic complex Lanczos algorithm                 //
//                  (without re-orthogonalization at each step)               //
//                 and storing each iteration information on disk             //
//                                                                            //
//                        last modification : 28/11/2003                      //
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


#include "config.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "GeneralTools/Endian.h"

#include <stdlib.h>
#include <string.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// directory  = directory where informations about Lanczos iterations have to be stored
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration

ComplexBasicLanczosAlgorithmWithDiskStorage::ComplexBasicLanczosAlgorithmWithDiskStorage(AbstractArchitecture* architecture, int nbrEigenvalue, 
											 int maxIter) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->NbrEigenvalue = nbrEigenvalue;
  this->V1 = ComplexVector();
  this->V2 = ComplexVector();
  this->V3 = ComplexVector();
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
  this->V1FileName = new char [8];
  this->V2FileName = new char [8];
  this->V3FileName = new char [8];
  strcpy(this->V1FileName, "vector1");
  strcpy(this->V2FileName, "vector2");
  strcpy(this->V3FileName, "vector3");
  this->LanczosFileName = new char [12];
  strcpy(this->LanczosFileName, "lanczos.dat");
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ComplexBasicLanczosAlgorithmWithDiskStorage::ComplexBasicLanczosAlgorithmWithDiskStorage(const ComplexBasicLanczosAlgorithmWithDiskStorage& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->V1FileName = new char [8];
  this->V2FileName = new char [8];
  this->V3FileName = new char [8];
  this->LanczosFileName = new char [12];
  strcpy(this->V1FileName, algorithm.V1FileName);
  strcpy(this->V2FileName, algorithm.V2FileName);
  strcpy(this->V3FileName, algorithm.V3FileName);
  strcpy(this->LanczosFileName, algorithm.LanczosFileName);
}

// destructor
//

ComplexBasicLanczosAlgorithmWithDiskStorage::~ComplexBasicLanczosAlgorithmWithDiskStorage() 
{
  delete[] this->V1FileName;
  delete[] this->V2FileName;
  delete[] this->V3FileName;
  delete[] this->LanczosFileName;  
}

// initialize Lanczos algorithm with a random vector
//

void ComplexBasicLanczosAlgorithmWithDiskStorage::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = ComplexVector (Dimension);
  this->V2 = ComplexVector (Dimension);
  this->V3 = ComplexVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    this->V1[i] = (rand() - 32767) * 0.5;
  this->V1 /= this->V1.Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ComplexBasicLanczosAlgorithmWithDiskStorage::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = vector;
  this->V2 = ComplexVector (Dimension);
  this->V3 = ComplexVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// resume Lanczos algorithm from disk datas in current directory
//

void ComplexBasicLanczosAlgorithmWithDiskStorage::ResumeLanczosAlgorithm()
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = ComplexVector (Dimension);
  this->V2 = ComplexVector (Dimension);
  this->V3 = ComplexVector (Dimension);
  this->ReadState();
}
  
// get last produced vector
//
// return value = reference on last produced vector

Vector& ComplexBasicLanczosAlgorithmWithDiskStorage::GetGroundState()
{
  return this->V2;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ComplexBasicLanczosAlgorithmWithDiskStorage::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V1, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->V1 * this->V2).Re;
      this->V2.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index), 
				    this->V1);
      this->V2 /= this->V2.Norm(); 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &this->V2, &this->V3);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
      this->V1.WriteVector(this->V1FileName);
      this->V2.WriteVector(this->V2FileName);
      this->V3.WriteVector(this->V3FileName);
      this->WriteState();
    }
  else
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  ComplexVector* TmpVector = new ComplexVector[2];
  double* TmpCoefficient = new double[2];
  //cout << this->Index << " " << Dimension << endl;
  for (int i = this->Index + 2; i < Dimension; ++i)
    {
      TmpVector[0] = this->V1;
      TmpVector[1] = this->V2;
      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
      AddComplexLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);
      this->V3 /= this->V3.Norm();
      this->V3.WriteVector(this->V3FileName);
      ComplexVector TmpV (this->V1);
      this->V1 = this->V2;
      this->V2 = this->V3;
      this->V3 = TmpV;
      char* TmpName = this->V1FileName;
      this->V1FileName = this->V2FileName;
      this->V2FileName = this->V3FileName;
      this->V3FileName = TmpName;
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation3 (this->Hamiltonian, &this->V2, &this->V3);
      Operation3.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = (this->V1 * this->V3).Re;
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = (this->V2 * this->V3).Re;
      this->V3.WriteVector(this->V3FileName);
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

bool ComplexBasicLanczosAlgorithmWithDiskStorage::TestConvergence ()
{
  if (((fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
	     (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1)))) ||
       (fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
	(this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1))))))
    return true;
  else
    return false;
}

// write current Lanczos state on disk
//
// return value = true if no error occurs

bool ComplexBasicLanczosAlgorithmWithDiskStorage::WriteState()
{
  ofstream File;
  File.open(this->LanczosFileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Index);
  WriteLittleEndian(File, this->PreviousLastWantedEigenvalue);
  WriteLittleEndian(File, this->EigenvaluePrecision);
  int TmpDimension = this->TridiagonalizedMatrix.GetNbrRow();
  WriteLittleEndian(File, TmpDimension);
  File.write(this->V1FileName, strlen(this->V1FileName) + 1);
  File.write(this->V2FileName, strlen(this->V2FileName) + 1);
  File.write(this->V3FileName, strlen(this->V3FileName) + 1);
  for (int i = 0; i <= (this->Index + 1); ++i)    
    {    
      WriteLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  for (int i = 0; i <= this->Index; ++i)
    {
      WriteLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  File.close();  
  return true;
}

// read current Lanczos state from disk
//
// return value = true if no error occurs

bool ComplexBasicLanczosAlgorithmWithDiskStorage::ReadState()
{
  ifstream File;
  File.open(this->LanczosFileName, ios::binary | ios::in);
  ReadLittleEndian(File, this->Index);
  ReadLittleEndian(File, this->PreviousLastWantedEigenvalue);
  ReadLittleEndian(File, this->EigenvaluePrecision);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
  this->V1FileName = new char [8];
  File.read(this->V1FileName, 8);
  this->V2FileName = new char [8];
  File.read(this->V2FileName, 8);
  this->V3FileName = new char [8];
  File.read(this->V3FileName, 8);
  this->V1.ReadVector(this->V1FileName);
  this->V2.ReadVector(this->V2FileName);
  this->V3.ReadVector(this->V3FileName);
  for (int i = 0; i <= (this->Index + 1); ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.DiagonalElement(i));
    }
  for (int i = 0; i <= this->Index; ++i)
    {
      ReadLittleEndian(File, this->TridiagonalizedMatrix.UpperDiagonalElement(i));
    }
  File.close();  
  this->Diagonalize();
  this->DiagonalizedMatrix.SortMatrixUpOrder();
  return true;
}


