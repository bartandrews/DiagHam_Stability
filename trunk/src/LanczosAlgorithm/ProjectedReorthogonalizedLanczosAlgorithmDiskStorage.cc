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


#include "LanczosAlgorithm/ProjectedReorthogonalizedLanczosAlgorithmDiskStorage.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"

#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <iostream>


using std::cout;
using std::endl;
using std::ios;


// flag switching whether to deallocate one vector upon application of the Hamiltonian in the main lanczos iteration
//#define DEALLOCATE_AT_H



// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ProjectedReorthogonalizedLanczosAlgorithmDiskStorage(AbstractHamiltonian** projectors, int nbrProjectors, AbstractArchitecture* architecture, 
													   int nbrEigenvalue, int maxIter, int nbrStorageVectors, int projectorIterMax,
													   double projectorPrecision, bool restartProjection, bool strongConvergence)
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MainIterMax = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->V1 = RealVector();
  this->V2 = RealVector();
  this->V3 = RealVector();

  this->VectorDimension = 0;
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MainIterMax, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MainIterMax, true);
    }
  else
    {
      this->MainIterMax = 1<<12;
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
  // projector fields:
  this->Projectors=projectors;
  this->NbrProjectors=nbrProjectors;
  this->InternalIndex=0;
  if (projectorIterMax > 0)
    {
      this->ProjectorIterMax=projectorIterMax;
      this->InternalTridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(ProjectorIterMax, true);
      this->InternalDiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(ProjectorIterMax);
    }
  else
    {
      this->ProjectorIterMax=1<<8;
      this->InternalTridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->InternalDiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->ProjectorPrecision=MACHINE_PRECISION;
  if (projectorPrecision != 0.0)
    this->ProjectorPrecision = projectorPrecision;
  this->RestartProjection = restartProjection;
  this->PreviousProjectorGroundstate=0.0;
  
  // initialize storage
  MainLanczosVectorFlags = new int[MainIterMax];
  for (int i=0; i< MainIterMax; ++i)
    MainLanczosVectorFlags[i]=0;
  ProjectorLanczosVectorFlags = new int[ProjectorIterMax];
  for (int i=0; i< ProjectorIterMax; ++i)
    ProjectorLanczosVectorFlags[i]=0;
  this->NbrStorageVectors = nbrStorageVectors;
  this->LanczosVectorStorage = new RealVector[NbrStorageVectors];

  this->VectorStorageFlags = new int[NbrStorageVectors];
  for (int i=0; i< NbrStorageVectors; ++i)
    VectorStorageFlags[i]=Empty;
  this->NbrStoredMain=0;
  this->NbrStoredInternal=0;
  this->Swap1Index = -1;
  this->Swap2Index = -1;
  this->VectorMarkers = BitMarkers(this->MainIterMax);
  this->TmpOutputName = new char[50];
  // various temporary arrays for RunLanczosAlgorithm
  this->OrthogonalizationSet = new RealVector*[2+NbrStorageVectors];
  this->OrthogonalizationCoef = new double[2+NbrStorageVectors];
  this->TmpCoefficient = new double[2];
  this->TmpScalarProduct = new double[2];
  this->TmpVectorArray = new RealVector[2];
  this->TmpVectorPtrArray = new RealVector*[2];
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ProjectedReorthogonalizedLanczosAlgorithmDiskStorage(const ProjectedReorthogonalizedLanczosAlgorithmDiskStorage& algorithm) 
{
  this->Index = algorithm.Index;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->V1 = algorithm.V1;
  this->V2 = algorithm.V2;
  this->V3 = algorithm.V3;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;
  this->VectorDimension = algorithm.VectorDimension;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->Projectors=algorithm.Projectors;
  this->NbrProjectors=algorithm.NbrProjectors;
  this->MainLanczosVectorFlags = algorithm.MainLanczosVectorFlags;
  this->MainIterMax = algorithm.MainIterMax;
  this->ProjectorLanczosVectorFlags = algorithm.ProjectorLanczosVectorFlags;
  this->ProjectorIterMax = algorithm.ProjectorIterMax;
  this->LanczosVectorStorage = algorithm.LanczosVectorStorage;
  this->NbrStorageVectors = algorithm.NbrStorageVectors;
  this->VectorStorageFlags = algorithm.VectorStorageFlags;
  this->NbrStoredMain = algorithm.NbrStoredMain;
  this->NbrStoredInternal = algorithm.NbrStoredInternal;
  this->VectorMarkers = algorithm.VectorMarkers;
  this->TmpOutputName = new char[50];
  this->OrthogonalizationSet = new RealVector*[2+NbrStorageVectors];
  this->OrthogonalizationCoef = new double[2+NbrStorageVectors];
  this->TmpCoefficient = new double[2];
  this->TmpScalarProduct = new double[2];
  this->TmpVectorArray = new RealVector[2];
  this->TmpVectorPtrArray = new RealVector*[2];
}

// destructor
//

ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::~ProjectedReorthogonalizedLanczosAlgorithmDiskStorage() 
{
  // clean up core of Projector
  this->ClearProjectorLanczosAlgorithm();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] this->VectorStorageFlags;
      delete [] this->MainLanczosVectorFlags;
      delete [] this->ProjectorLanczosVectorFlags;
      if (NbrStorageVectors>0)
	delete [] this->LanczosVectorStorage;
    }
  delete [] this->PreviousWantedEigenvalues;
  delete [] this->OrthogonalizationSet;
  delete [] this->OrthogonalizationCoef;  
  delete [] this->TmpCoefficient;
  delete [] this->TmpScalarProduct;
  delete [] this->TmpVectorArray;
  delete [] this->TmpVectorPtrArray;
  delete [] TmpOutputName;
}

// initialize Lanczos algorithm with a random vector
//

void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::InitializeLanczosAlgorithm() 
{
  if (this->Hamiltonian == NULL)
    {
      cout << "A Hamiltonian has to be declared before initializing a Lanczos Algorithm"<<endl;
      exit(-1);
    }
  this->VectorDimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->V1 = RealVector (VectorDimension);
  this->V2 = RealVector (VectorDimension);
  this->V3 = RealVector (VectorDimension);
  int Shift = RAND_MAX / 2;
  double Scale = 1.0 / ((double) Shift);
  for (int i = 0; i < VectorDimension; i++)
    {
      this->V1[i] = Scale * ((double) (rand() - Shift));
    }
  this->V1 /= this->V1.Norm();
  
  // project initial vector
  for (int p=0; p<NbrProjectors; ++p)
    this->ProjectVector(p);
  //cout << "calling this->SaveVector(V1,0) on line "<<__LINE__<<endl;
  this->SaveVector(this->V1,0,true,true);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
  for (int i=0; i<NbrStorageVectors; ++i)
    this->LanczosVectorStorage[i].Resize(VectorDimension);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::InitializeLanczosAlgorithm(const Vector& vector) 
{
  if (this->Hamiltonian == NULL)
    {
      cout << "A Hamiltonian has to be declared before initializing a Lanczos Algorithm"<<endl;
      exit(-1);
    }
  this->VectorDimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if (VectorDimension != vector.GetVectorDimension())
    {
      cout << "initial vector does not match dimension of Hilbert-space"<<endl;
    }
  this->V1 = vector;
  this->V2 = RealVector (VectorDimension);
  this->V3 = RealVector (VectorDimension);
  // project initial vector
  for (int p=0; p<NbrProjectors; ++p)
    this->ProjectVector(p);
  //cout << "calling this->SaveVector(V1,0) on line "<<__LINE__<<endl;
  this->SaveVector(this->V1,0,true,true);

  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
  for (int i=0; i<NbrStorageVectors; ++i)
    this->LanczosVectorStorage[i].Resize(VectorDimension);
}

// resume Lanczos algorithm from disk datas in current directory
//

void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ResumeLanczosAlgorithm()
{
  this->ReadState();
  this->VectorDimension = this->V1.GetVectorDimension();
  if (this->Hamiltonian != NULL)
    {
      if (this->Hamiltonian->GetHilbertSpaceDimension() != this->VectorDimension)
	{
	  cout << "Hamiltonian does not match dimension of stored vectors in resuming"<<endl;
	  exit(-1);
	}
    }
  for (int i=0; i<NbrStorageVectors; ++i)
    this->LanczosVectorStorage[i].Resize(VectorDimension);
  for (int i=0; i<this->Index + 2; ++i)
    MainLanczosVectorFlags[i]|=SavedOnDisk;
}
  
// get last produced vector
//
// return value = reference on lest produced vector

Vector& ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::GetGroundState()
{
  Vector *TmpVector = this->GetEigenstates(1);
  return *TmpVector;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::GetEigenstates(int nbrEigenstates)
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
  RealVector **TmpVectors = this->OrthogonalizationSet;
  int CollectionSize=3;
  TmpVectors[0]=&(this->V1);
  TmpVectors[1]=&(this->V2);
  TmpVectors[2]=&(this->V3);
  // initialize eigenstates
  for (int i = 0; i < nbrEigenstates; ++i)
    Eigenstates[i] = RealVector (this->VectorDimension, true);
  int CurrentNbrVector=0;
  int StartIndex=0;
  int CollectionIndex=0;
  
  while (CurrentNbrVector<this->TridiagonalizedMatrix.GetNbrRow())
    {
      this->ReadVector(*(TmpVectors[CollectionIndex]),CurrentNbrVector);
      ++CollectionIndex;
      ++CurrentNbrVector;

      if (CollectionIndex==CollectionSize)
	{
	  for (int i = 0; i < nbrEigenstates; ++i)
	    {
	      CollectionIndex=0;
	      for (int k=StartIndex; k<CurrentNbrVector; ++k, ++CollectionIndex)
		TmpCoefficents[CollectionIndex] = TmpEigenvector(k, i);
	      AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), TmpVectors, CollectionSize, TmpCoefficents);
	      Operation.ApplyOperation(this->Architecture);
	    }
	  StartIndex=CurrentNbrVector;
	  CollectionIndex=0;
	}
    }
  if (CollectionIndex>0)
    {
      for (int i = 0; i < nbrEigenstates; ++i)
	{
	  CollectionIndex=0;
	  for (int k=StartIndex; k<CurrentNbrVector; ++k, ++CollectionIndex)
	    TmpCoefficents[CollectionIndex] = TmpEigenvector(k, i);
	  AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), TmpVectors, CollectionIndex, TmpCoefficents);
	  Operation.ApplyOperation(this->Architecture);
	}
    }
  for (int i = 0; i < nbrEigenstates; ++i)
    Eigenstates[i] /= Eigenstates[i].Norm();

  delete[] TmpCoefficents;
  return Eigenstates;
}



// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->TridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->TridiagonalizedMatrix.GetNbrRow() + 2;
      this->TridiagonalizedMatrix.Resize(Dimension, Dimension);
      // begin step 0->1
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V1, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.DiagonalElement(Index) = (this->V1 * this->V2);
      this->V2.AddLinearCombination(-this->TridiagonalizedMatrix.DiagonalElement(this->Index), this->V1);
      this->V2 /= this->V2.Norm();
      
      // apply projector here
      // (V1 already saved in initialization)
      // swap V1, V2
      {
	RealVector TmpV (this->V1);
	this->V1 = this->V2;
	this->V2 = TmpV; 
      }
      //cout << "Line "<<__LINE__<<": V1="<<endl<<V1;
      bool RequireReorthogonalization = false;
      for (int p=0; p<NbrProjectors; ++p)
	RequireReorthogonalization |= this->ProjectVector(p);
      // swap V1, V2
      {
	RealVector TmpV (this->V1);
	this->V1 = this->V2;
	this->V2 = TmpV; 
      }
      //cout << "calling this->ReadVector(V1,0) on line "<<__LINE__<<endl;
      this->ReadVector(V1,0);      
      
      // end projector

      // orthogonalizing again, if needed
      if (RequireReorthogonalization)
	{
	  double NewScalarProd = this->V1 * this->V2;
	  this->V2.AddLinearCombination(-NewScalarProd, this->V1);
	  this->V2 /= this->V2.Norm();
	}
      //cout << "calling this->SaveVector(V2,1) on line "<<__LINE__<<endl;
      this->SaveVector(V2,1,true,true);
      // ============ end step 0->1 ===========
      // begin step 1->2
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
  double TmpScalarProduct[2];
  for (int i = this->Index + 2; i < Dimension; i++)
    {
      TmpVectorArray[0] = this->V1;
      TmpVectorArray[1] = this->V2;
      TmpCoefficient[0] = -this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index);
      TmpCoefficient[1] = -this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1);
      AddRealLinearCombinationOperation Operation4 (&(this->V3),  TmpVectorArray, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);	  
      this->V3 /= this->V3.Norm();
      
      // apply projector here
      // vector0..i-1 are already saved
      // swap V1, V3
      {
	RealVector TmpV (this->V1);
	this->V1 = this->V3;
	this->V3 = TmpV;
      }
      bool RequireReorthogonalization = false;
      for (int p=0; p<NbrProjectors; ++p)
	RequireReorthogonalization |= this->ProjectVector(p);
      // full reorthogonalization procedure on vectors 0 ... i-3
      this->FullReorthogonalisation(i-2);
      
      // swap V1, V3
      {
	RealVector TmpV (this->V1);
	this->V1 = this->V3;
	this->V3 = TmpV; 
      }
      //cout << "calling this->ReadVector(V1,"<<i-2<<") on line "<<__LINE__<< " flags="<<std::hex<<MainLanczosVectorFlags[i-2]<<std::dec<<endl;
      this->ReadVector(V1,i-2);
      //cout << "calling this->ReadVector(V2,"<<i-1<<") on line "<<__LINE__<<" flags="<<std::hex<<MainLanczosVectorFlags[i-1]<<std::dec<<endl;
      this->ReadVector(V2,i-1);
      // end projector
	  
      // reorthogonalize once more, if needed:
      if (RequireReorthogonalization)
	{
	  // recalculate scalar products
	  RealVector* TmpVectorPtrArray[2];
	  TmpVectorPtrArray[0] = &(this->V1);
	  TmpVectorPtrArray[1] = &(this->V2);
	  MultipleRealScalarProductOperation Operation (&(this->V3), TmpVectorPtrArray, 2, TmpScalarProduct);
	  Operation.ApplyOperation(this->Architecture);
	  // perform subtractions
	  TmpVectorArray[0] = this->V1;
	  TmpVectorArray[1] = this->V2;
	  TmpCoefficient[0] = -TmpScalarProduct[0];
	  TmpCoefficient[1] = -TmpScalarProduct[1];
	  AddRealLinearCombinationOperation Operation1 (&(this->V3),  TmpVectorArray, 2, TmpCoefficient);
	  Operation1.ApplyOperation(this->Architecture);	  
	  this->V3 /= this->V3.Norm();
	}
      //cout << "calling this->SaveVector(V3,"<<i<<") on line "<<__LINE__<<endl;
      this->SaveVector(V3, i, true, true);
      this->WriteState();

#ifdef DEALLOCATE_AT_H
      RealVector TmpV (this->V2);
      this->V2 = this->V3;
      this->V3 = TmpV;	  
      this->V1 = RealVector();
#else
      RealVector TmpV (this->V1);
      this->V1 = this->V2;
      this->V2 = this->V3;
      this->V3 = TmpV;
#endif
      this->Index++;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &this->V2, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
#ifdef DEALLOCATE_AT_H
      //cout << "calling this->ReadVector(V1,"<<i-1<<") on line "<<__LINE__<<endl;
      this->ReadVector(V1,i-1);
#endif
      RealVector* TmpVectorPtrArray[2];
      TmpVectorPtrArray[0] = &(this->V1);
      TmpVectorPtrArray[1] = &(this->V2);
      MultipleRealScalarProductOperation Operation2 (&(this->V3), TmpVectorPtrArray, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->TridiagonalizedMatrix.UpperDiagonalElement(this->Index) = TmpScalarProduct[0];
      this->TridiagonalizedMatrix.DiagonalElement(this->Index + 1) = TmpScalarProduct[1];
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
  this->WriteState();
}


  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::TestConvergence ()
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

bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::WriteState()
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

bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ReadState()
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
  sprintf(TmpOutputName, "vector.%d", this->Index);
  this->V1.ReadVector(TmpOutputName);
  sprintf(TmpOutputName, "vector.%d", (this->Index + 1));
  this->V2.ReadVector(TmpOutputName);
  sprintf(TmpOutputName, "vector.%d", (this->Index + 2));
  this->V3.ReadVector(TmpOutputName);
  this->Diagonalize();
  this->DiagonalizedMatrix.SortMatrixUpOrder();
  return true;
}

// proceed through reorthogonalization of vectors up to given step
// nbrStep = number of last vector to reorthogonalize with
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::FullReorthogonalisation(int nbrStep)
{
  if (nbrStep<0) return;
  // first, orthogonalize with respect to vectors that are still held in memory
  int VectorFlags;
  int CollectionSize=0;
  VectorMarkers.ResetMarked();
  for (int i=0; i<nbrStep; ++i)
    {
      VectorFlags = MainLanczosVectorFlags[i];
      if (VectorFlags & SavedInMemory)
	{
	  //cout << "Orthogonalizing from memory with state "<<i<<endl;
	  VectorMarkers.MarkBit(i);
	  OrthogonalizationSet[CollectionSize]=&(LanczosVectorStorage[VectorFlags & VectorIndexMask]);
	  ++CollectionSize;
	}
    }
  if (CollectionSize>0) // found any vectors in memory?
    {
      MultipleRealScalarProductOperation Operation (&(this->V1), OrthogonalizationSet, CollectionSize, OrthogonalizationCoef);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < CollectionSize; j++)
	OrthogonalizationCoef[j] *= -1.0;
      AddRealLinearCombinationOperation Operation2 (&(this->V1), OrthogonalizationSet, CollectionSize, OrthogonalizationCoef);
      Operation2.ApplyOperation(this->Architecture);
    }
  // treat remaining vectors
  // define storage space
  OrthogonalizationSet[0] = &(this->V2);
  OrthogonalizationSet[1] = &(this->V3);
  CollectionSize=2;
  // claim additional storage from any storage of internal vectors
  if (NbrStorageVectors-NbrStoredMain>0)
    {
      //cout << "May claim " << NbrStorageVectors-NbrStoredMain << " empty/internal vectors!"<<endl;
      for (int i=0; i<NbrStorageVectors; ++i)
	{
	  if (VectorStorageFlags[i] & ProjectorLanczosVector)
	    {
	      //cout << "claiming internal vector: "<<i<<" with dimension "<<LanczosVectorStorage[i].GetVectorDimension()<<endl;
	      OrthogonalizationSet[CollectionSize] = &(LanczosVectorStorage[i]);
	      ++CollectionSize;	      
	      ProjectorLanczosVectorFlags[VectorStorageFlags[i] & VectorIndexMask] = 0;
	      --NbrStoredInternal;
	      VectorStorageFlags[i]=Empty;
	    }
	  else if (VectorStorageFlags[i]==Empty)
	    {
	      //cout << "claiming empty vector: "<<i<<" with dimension "<<LanczosVectorStorage[i].GetVectorDimension()<<endl;
	      OrthogonalizationSet[CollectionSize] = &(LanczosVectorStorage[i]);
	      ++CollectionSize;
	    }
	}
    }
  int CurrentNbrVector=0;
  int CollectionIndex=0;
  while (CurrentNbrVector<nbrStep)
    {
      if (VectorMarkers.IsNotMarked(CurrentNbrVector))
	{
	  sprintf(TmpOutputName, "vector.%d", CurrentNbrVector);
	  OrthogonalizationSet[CollectionIndex++]->ReadVector(TmpOutputName);
	}
      if (CollectionIndex==CollectionSize)
	{
	  MultipleRealScalarProductOperation Operation (&(this->V1), OrthogonalizationSet, CollectionSize, OrthogonalizationCoef);
	  Operation.ApplyOperation(this->Architecture);
	  for (int j = 0; j < CollectionSize; j++)
	    OrthogonalizationCoef[j] *= -1.0;
	  AddRealLinearCombinationOperation Operation2 (&(this->V1), OrthogonalizationSet, CollectionSize, OrthogonalizationCoef);
	  Operation2.ApplyOperation(this->Architecture);
	  CollectionIndex=0;
	}
      ++CurrentNbrVector;
    }
  if (CollectionIndex>0)
    {
      MultipleRealScalarProductOperation Operation (&(this->V1), OrthogonalizationSet, CollectionIndex, OrthogonalizationCoef);
      Operation.ApplyOperation(this->Architecture);
      for (int j = 0; j < CollectionSize; j++)
	OrthogonalizationCoef[j] *= -1.0;
      AddRealLinearCombinationOperation Operation2 (&(this->V1), OrthogonalizationSet, CollectionIndex, OrthogonalizationCoef);
      Operation2.ApplyOperation(this->Architecture);
      CollectionIndex=0;
    }
}



// project the vector stored in V1 to the groundstate of the given Projector number
// nbrProjector = number of projector to use for this run
// return = true if Lanczos was run, or false if no projection needed to be done
bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ProjectVector(int nbrProjector)
{
  this->InitializeProjectorLanczosAlgorithm();

  cout << "[ Projection "<<nbrProjector<<"...";
  bool Done = this->RunProjectorLanczosAlgorithm(nbrProjector, 3);
  if (Done)
    {
      cout << "not required! ]"<<endl;
      return false;
    }
  else
    {      
      int CurrentNbrIterLanczos = 4;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);
      while ((this->TestProjectorConvergence() == false) && (CurrentNbrIterLanczos < this->ProjectorIterMax))
	{
	  cout << endl;
	  ++CurrentNbrIterLanczos;
	  this->RunProjectorLanczosAlgorithm(nbrProjector, 1);
	  gettimeofday (&(TotalEndingTime), 0);
	  double Precision = fabs((PreviousProjectorGroundstate - ProjectorGroundstate) / PreviousProjectorGroundstate);
	  double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);		      
	  TotalStartingTime.tv_usec = TotalEndingTime.tv_usec;
	  TotalStartingTime.tv_sec = TotalEndingTime.tv_sec;
	  cout << "   "<<CurrentNbrIterLanczos<<": " << ProjectorGroundstate << " " << Precision << " ("<<Dt<<" s)";
	}
      cout << "... DONE ]"<<endl;
      bool Restart=false;
      if (this->RestartProjection)
	Restart = !(this->TestProjectorConvergence());
      this->GetProjectorGroundState();
      if (Restart)
	this->ProjectVector(nbrProjector);
      return true;
    }
}

// test if convergence has been reached in internal lanczos
//
// return value = true if convergence has been reached

bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::TestProjectorConvergence ()
{
  if ( (fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0)) < this->ProjectorPrecision) ||
       (fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0) - this->PreviousProjectorGroundstate) < 
	(this->ProjectorPrecision * fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0))) ) )
    return true;
  else
    return false;
}


// run current internal Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrProjector = number of projector to use for this run
// nbrIter = number of iteration to do 

bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::RunProjectorLanczosAlgorithm (int nbrProjector, int nbrIter) 
{
  int Dimension;
  if (this->InternalIndex == 0)
    {
      Dimension = this->InternalTridiagonalizedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->InternalTridiagonalizedMatrix.GetNbrRow() + 2;
      this->InternalTridiagonalizedMatrix.Resize(Dimension, Dimension);
      VectorHamiltonianMultiplyOperation Operation1 (this->Projectors[nbrProjector], &this->V1, &this->V2);
      Operation1.ApplyOperation(this->Architecture);
      this->InternalTridiagonalizedMatrix.DiagonalElement(InternalIndex) = (this->V1 * this->V2);
      if (fabs(this->InternalTridiagonalizedMatrix.DiagonalElement(InternalIndex))<this->ProjectorPrecision)
	{
	  return true;
	}
      // running internal Lanczos -> save starting vector, now
      //cout << "calling this->SaveVector(V1,0, interal) on line "<<__LINE__<<endl;
      this->SaveVector(V1, 0, false, true); // need to keep V1 in process
      this->V2.AddLinearCombination(-this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex), 
				    this->V1);
      this->V2 /= this->V2.Norm();
      //cout << "calling this->SaveVector(V2,"<<1<<",false,true) on line "<<__LINE__<<endl;
      this->SaveVector(V2, 1, false, true); // need to keep V2 in process, here
      // it may be necessary to resize V3 if this vector has been dropped earlier in main lanczos
      if (V3.GetVectorDimension()==0)
	V3.Resize(V2.GetVectorDimension());
      VectorHamiltonianMultiplyOperation Operation2 (this->Projectors[nbrProjector], &this->V2, &this->V3);
      Operation2.ApplyOperation(this->Architecture);
      this->InternalTridiagonalizedMatrix.UpperDiagonalElement(this->InternalIndex) = (this->V1 * this->V3);
      this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex + 1) = (this->V2 * this->V3);
    }
  else
    {
      Dimension = this->InternalTridiagonalizedMatrix.GetNbrRow() + nbrIter;
      this->InternalTridiagonalizedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->InternalIndex + 2; i < Dimension; i++)
    {
      TmpVectorArray[0] = this->V1;
      TmpVectorArray[1] = this->V2;
      TmpCoefficient[0] = -this->InternalTridiagonalizedMatrix.UpperDiagonalElement(this->InternalIndex);
      TmpCoefficient[1] = -this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex + 1);
      AddRealLinearCombinationOperation Operation4 (&(this->V3),  TmpVectorArray, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);
      this->V3 /= this->V3.Norm();
      //cout << "calling this->SaveVector(V3,"<<i<<",false) on line "<<__LINE__<<endl;
      this->SaveVector(V3, i, false, true); // need V3 in process
      
      if (false) // hardwired option "fast-disk" which forgets one of the vectors at this point to save memory
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
      ++this->InternalIndex;
      VectorHamiltonianMultiplyOperation Operation1 (this->Projectors[nbrProjector], &this->V2, &this->V3);
      Operation1.ApplyOperation(this->Architecture);
      // if (this->DiskFlag == true)
      // again, hardwired
      if (false)
	{
	  this->ReadVector(V1, i-1, false); // reread vector, and keep a copy
	}
      TmpVectorPtrArray[0] = &(this->V1);
      TmpVectorPtrArray[1] = &(this->V2);
      MultipleRealScalarProductOperation Operation2 (&(this->V3), TmpVectorPtrArray, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->InternalTridiagonalizedMatrix.UpperDiagonalElement(this->InternalIndex) = TmpScalarProduct[0];
      this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex + 1) = TmpScalarProduct[1];
    }
  if (this->PreviousProjectorGroundstate != 0.0)
    {
      this->PreviousProjectorGroundstate = this->InternalDiagonalizedMatrix.DiagonalElement(0);
      this->ProjectorDiagonalize();
      this->InternalDiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->ProjectorDiagonalize();
      this->InternalDiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousProjectorGroundstate = 2.0 * this->InternalDiagonalizedMatrix.DiagonalElement(0);
    }
  return false;
}

// get last vector produced by projector lanczos; result is stored in internal vector V1
//
// return value = reference on produced vector

Vector& ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::GetProjectorGroundState()
{
  RealMatrix TmpEigenvector (this->InternalTridiagonalizedMatrix.GetNbrRow(), this->InternalTridiagonalizedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->InternalTridiagonalizedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;
  
  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->InternalTridiagonalizedMatrix.GetNbrRow());
  SortedDiagonalizedMatrix.Copy(this->InternalTridiagonalizedMatrix);
  SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
  double* TmpComponents = new double [this->InternalTridiagonalizedMatrix.GetNbrRow()];
  for (int j = 0; j < this->InternalTridiagonalizedMatrix.GetNbrRow(); ++j)
    TmpComponents[j] = TmpEigenvector(j, 0);

  this->ReadVector(this->V2, this->InternalDiagonalizedMatrix.GetNbrRow()-1, false, false);
  this->V1.Copy(this->V2, TmpComponents[this->InternalDiagonalizedMatrix.GetNbrRow()-1]);
  for (int i = this->InternalDiagonalizedMatrix.GetNbrRow()-2; i > -1; --i)
    {
      this->ReadVector(V2, i, false, false);    
      this->V1.AddLinearCombination(TmpComponents[i], this->V2);
      cout << i << "/" << this->InternalDiagonalizedMatrix.GetNbrRow() << "           \r";
      cout.flush();
    }
  cout << endl;
  this->V1 /= this->V1.Norm();  
  delete[] TmpComponents;
  return this->V1;
}

// initialize Lanczos algorithm starting from the vector currently held in V1
//
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::InitializeProjectorLanczosAlgorithm()
{
  if (this->InternalIndex>0)
    this->ClearProjectorLanczosAlgorithm();
  // normalization and external orthogonalization should already have been done before calling this...
  // this->V1 /= this->V1.Norm();  
  this->InternalTridiagonalizedMatrix.ResizeAndClean(0, 0);
}


// clear storage associated with internal Lanczos algorithm
//
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ClearProjectorLanczosAlgorithm()
{
  // forget all previously stored internal lanczos vectors
  for (int i=0; i<NbrStorageVectors; ++i)
    if(VectorStorageFlags[i]&ProjectorLanczosVector)
      {
	VectorStorageFlags[i]=Empty;
	ProjectorLanczosVectorFlags[VectorStorageFlags[i]&StorageIndexMask]=0;
      }
  for (int i=0; i<InternalIndex+2; ++i)
    {
      //       cout << "ProjectorLanczosVectorFlags["<<i<<"]="<<hex<<ProjectorLanczosVectorFlags[i]<<dec<<endl;
      //       cout.setf(ios::dec, ios::basefield); // alternative way of setting output format
      if (ProjectorLanczosVectorFlags[i]&SavedOnDisk)
	{
	  sprintf(TmpOutputName,"vector-int.%d",i);
	  if( remove(TmpOutputName) != 0 )
	    cout<<"Error deleting file "<<TmpOutputName<<endl;
	  // 	  else
	  // 	    cout<<"Deleted file "<<TmpOutputName<<endl;
	}
      ProjectorLanczosVectorFlags[i]=0;
    }
  this->NbrStoredInternal=0;
  this->InternalIndex=0;
  this->PreviousProjectorGroundstate=0.0;
}

// diagonalize internal tridiagonalized matrix and find ground state energy
//
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ProjectorDiagonalize ()
{
  int Dimension = this->InternalTridiagonalizedMatrix.GetNbrRow();
  this->InternalDiagonalizedMatrix.Copy(this->InternalTridiagonalizedMatrix);
  this->InternalDiagonalizedMatrix.Diagonalize(50);
  this->ProjectorGroundstate = this->InternalDiagonalizedMatrix.DiagonalElement(0);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (this->InternalDiagonalizedMatrix.DiagonalElement(DiagPos) < this->ProjectorGroundstate)
      this->ProjectorGroundstate = this->InternalDiagonalizedMatrix.DiagonalElement(DiagPos);  
  return;
}


// save vector - test version, always saving to disk - SLOWER!
// vec = vector to be savec
// index = vector index in Lanczos routine
// mainLanczos = flag indicating whether vector is part of main lanczos algorithm
// keepOriginal = flag indicating whether original vector needs to be kept in place
// return = true on success
bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::SaveVectorTest(RealVector &vec, int index, bool mainLanczos, bool keepOriginal)
{
  int StorageFlag=index & StorageIndexMask;
  int VectorFlag=0;

  if (mainLanczos==true)
    {
      StorageFlag |= MainLanczosVector;
      sprintf(TmpOutputName,"vector.%d",index);
      if (MainLanczosVectorFlags[index] & SavedOnDisk)
	{
	  cout << "attention, overwriting vector "<<TmpOutputName<<endl;
	}
      vec.WriteVector(TmpOutputName);
      VectorFlag |= SavedOnDisk;
      MainLanczosVectorFlags[index] = VectorFlag;
    }
  else
    {
      StorageFlag |= ProjectorLanczosVector;
      sprintf(TmpOutputName,"vector-int.%d",index);
      if (ProjectorLanczosVectorFlags[index] & SavedOnDisk)
	{
	  cout << "attention, overwriting vector "<<TmpOutputName<<endl;
	}
      vec.WriteVector(TmpOutputName);
      VectorFlag |= SavedOnDisk;
      ProjectorLanczosVectorFlags[index] = VectorFlag;
    }
  return true;
}


// reread vector - test version, always reading from disk
// vec = vector to be retrieved
// index = vector index in Lanczos routine
// mainLanczos = flag indicating whether vector is part of main lanczos algorithm
// keepCopy = flag indicating whether the saved vector still needs to be kept in memory after reloading
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ReadVectorTest(RealVector &vec, int index, bool mainLanczos, bool keepCopy)
{
  int VectorFlags;
  if (mainLanczos == true)
    VectorFlags = MainLanczosVectorFlags[index];
  else
    VectorFlags = ProjectorLanczosVectorFlags[index];
  if (VectorFlags==0)
    {
      cout << "Error: cannot retrieve vector "<<index<<" for "<<(mainLanczos?"main":"internal")<<" Lanczos!"<<endl;
      exit(-1);
    }
  if (VectorFlags & SavedOnDisk)
    {
      if (mainLanczos == true)
	sprintf(TmpOutputName, "vector.%d", index);
      else
	sprintf(TmpOutputName, "vector-int.%d", index);
      vec.ReadVector(TmpOutputName);
    }
  return;
}



// save vector, either to internal memory, or disk, or both
// vec = vector to be savec
// index = vector index in Lanczos routine
// mainLanczos = flag indicating whether vector is part of main lanczos algorithm
// keepOriginal = flag indicating whether original vector needs to be kept in place
// return = true on success
bool ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::SaveVector(RealVector &vec, int index, bool mainLanczos, bool keepOriginal)
{
  int StorageFlag=index & StorageIndexMask;
  int VectorFlag=0;
  
  if (mainLanczos==true)
    {
      StorageFlag |= MainLanczosVector;
      
      // always write to disk here
      sprintf(TmpOutputName,"vector.%d",index);
      if (MainLanczosVectorFlags[index] & SavedOnDisk)
	{
	  cout << "attention, overwriting vector "<<TmpOutputName<<endl;
	}
      vec.WriteVector(TmpOutputName);
      VectorFlag |= SavedOnDisk;
      MainLanczosVectorFlags[index] = VectorFlag;
      if (NbrStorageVectors>0)
	{
	  // test whether this vector is already held in storage
	  if (MainLanczosVectorFlags[index] & SavedInMemory)
	    {
	      for (int i=0; i<NbrStorageVectors; ++i)
		if ((VectorStorageFlags[i] & StorageIndexMask) == index) // find vector
		  {
		    if (keepOriginal)
		      {
			LanczosVectorStorage[i].Copy(vec);
			cout << "attention, overwriting vector "<<index<<" in RAM"<<endl;
		      }
		    else
		      {
			RealVector TmpV(vec);
			vec=LanczosVectorStorage[i];
			LanczosVectorStorage[i]=TmpV;
		      }
		    break;
		  }
	    }
	  else
	    {
	      // then search for additional storage in live memory
	      if (NbrStoredMain<NbrStorageVectors) // have more space to give away?
		{
		  for (int i=0; i<NbrStorageVectors; ++i)
		    {
		      if (VectorStorageFlags[i]==Empty) // new empty position?
			{
			  ++NbrStoredMain;
			  VectorStorageFlags[i] = StorageFlag;
			  if (keepOriginal)
			    {
			      LanczosVectorStorage[i].Copy(vec);
			    }
			  else
			    {
			      RealVector TmpV(vec);
			      vec=LanczosVectorStorage[i];
			      LanczosVectorStorage[i]=TmpV;
			    }
			  VectorFlag |= (i & VectorIndexMask);
			  VectorFlag |= SavedInMemory;
			  MainLanczosVectorFlags[index]=VectorFlag;
			  break;
			}
		      if ((VectorStorageFlags[i]&ProjectorLanczosVector)!=0) // internal storage that can be overwritten?
			{
			  ++NbrStoredMain;
			  --NbrStoredInternal;
			  // mark state as forgotten
			  int TmpIndex = (VectorStorageFlags[i]&StorageIndexMask);
			  ProjectorLanczosVectorFlags[TmpIndex]=0;
			  VectorStorageFlags[i] = StorageFlag;
			  if (keepOriginal)
			    {
			      LanczosVectorStorage[i].Copy(vec);
			    }
			  else
			    {
			      RealVector TmpV(vec);
			      vec=LanczosVectorStorage[i];
			      LanczosVectorStorage[i]=TmpV;
			    }
			  VectorFlag |= (i & VectorIndexMask);
			  VectorFlag |= SavedInMemory;
			  MainLanczosVectorFlags[index]=VectorFlag;
			  break;
			}
		    }
		}
	      else // have no more space to give away: overwrite oldest vector
		{
		  int MinPos=0;
		  int MinIndex=MainIterMax;
		  for (int i=0; i<NbrStorageVectors; ++i)
		    {
		      if ((VectorStorageFlags[i]&StorageIndexMask) < MinIndex)
			{
			  MinPos=i;
			  MinIndex=(VectorStorageFlags[i]&StorageIndexMask);
			}
		    }
		  // check for consistence: state should have been written to disk
		  if ((MainLanczosVectorFlags[MinIndex] & SavedOnDisk) == 0)
		    {
		      cout << "Problem in ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::SaveVector - should have saved Main vector no. MinIndex, but not flagged"<<endl;
		      exit(-1);
		    }
		  // mark state as forgotten from memory
		  MainLanczosVectorFlags[MinIndex]&= ~SavedInMemory;
		  MainLanczosVectorFlags[MinIndex]&= ~VectorIndexMask;
		  VectorStorageFlags[MinPos] = StorageFlag;
		  if (keepOriginal)
		    {
		      LanczosVectorStorage[MinPos].Copy(vec);
		    }
		  else
		    {
		      RealVector TmpV(vec);
		      vec=LanczosVectorStorage[MinPos];
		      LanczosVectorStorage[MinPos]=TmpV;
		    }
		  VectorFlag |= (MinPos & VectorIndexMask);
		  VectorFlag |= SavedInMemory;
		  MainLanczosVectorFlags[index]=VectorFlag;
		}
	    }
	}
    }
  else // store vector for internal lanczos
    {
      StorageFlag |= ProjectorLanczosVector;
      if (NbrStorageVectors>0) // have internal storage at all?
	{
	  // search for additional storage in live memory
	  if (NbrStoredMain+NbrStoredInternal<NbrStorageVectors) // have more space to give away?
	    {
	      for (int i=0; i<NbrStorageVectors; ++i)
		if (VectorStorageFlags[i]==Empty) // new empty position?
		  {
		    ++NbrStoredInternal;
		    VectorStorageFlags[i] = StorageFlag;
		    if (keepOriginal)
		      {
			LanczosVectorStorage[i].Copy(vec);
		      }
		    else
		      {
			RealVector TmpV(vec);
			vec=LanczosVectorStorage[i];
			LanczosVectorStorage[i]=TmpV;
		      }
		    VectorFlag |= (i & VectorIndexMask);
		    VectorFlag |= SavedInMemory;
		    ProjectorLanczosVectorFlags[index]=VectorFlag;
		    break;
		  }
	    }
	  else // have no more space to give away
	    {
	      if (NbrStoredMain>2) // have more than two last vectors from main Lanczos -> overwrite oldest
		{
		  int MinPos=0;
		  int MinIndex=MainIterMax;
		  for (int i=0; i<NbrStorageVectors; ++i)
		    {
		      if (((VectorStorageFlags[i]&MainLanczosVector) != 0)
			  &&((VectorStorageFlags[i]&StorageIndexMask) < MinIndex))
			{
			  MinPos=i;
			  MinIndex=(VectorStorageFlags[i]&StorageIndexMask);
			}
		    }
		  // check for consistence: state should have been written to disk
		  if ((MainLanczosVectorFlags[MinIndex] & SavedOnDisk) == 0)
		    {
		      cout << "Problem in ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::SaveVector - should have saved Main vector no. MinIndex, but not flagged"<<endl;
		      exit(-1);
		    }
		  // mark state as forgotten
		  MainLanczosVectorFlags[MinIndex] &= ~SavedInMemory;
		  MainLanczosVectorFlags[MinIndex] &= ~VectorIndexMask;
		  --NbrStoredMain;
		  VectorStorageFlags[MinPos] = StorageFlag;
		  if (keepOriginal)
		    {
		      LanczosVectorStorage[MinPos].Copy(vec);
		    }
		  else
		    {
		      RealVector TmpV(vec);
		      vec=LanczosVectorStorage[MinPos];
		      LanczosVectorStorage[MinPos]=TmpV;
		    }
		  VectorFlag |= (MinPos & VectorIndexMask);
		  VectorFlag |= SavedInMemory;
		  ProjectorLanczosVectorFlags[index]=VectorFlag;
		  ++NbrStoredInternal;
		}
	      else // have two or less states from main lanczos
		{
		  // did we store any internal vectors at all? -> write oldest internal state to disk
		  if (NbrStoredInternal>0)
		    {
		      int MinPos=0;
		      int MinIndex=ProjectorIterMax;
		      for (int i=0; i<NbrStorageVectors; ++i)
			{
			  if (((VectorStorageFlags[i]&ProjectorLanczosVector) != 0)
			      &&((VectorStorageFlags[i]&StorageIndexMask) < MinIndex))
			    {
			      MinPos=i;
			      MinIndex=(VectorStorageFlags[i]&StorageIndexMask);
			    }
			}
		      // save oldest vector to disk
		      sprintf(TmpOutputName,"vector-int.%d",MinIndex);
		      //cout << "Writing to disk: "<<TmpOutputName<<endl;
		      LanczosVectorStorage[MinPos].WriteVector(TmpOutputName);
		      ProjectorLanczosVectorFlags[MinIndex] |= SavedOnDisk;
		      ProjectorLanczosVectorFlags[MinIndex] &= ~SavedInMemory;
		      ProjectorLanczosVectorFlags[MinIndex] &= ~VectorIndexMask;

		      // save new vector in memory
		      VectorStorageFlags[MinPos] = StorageFlag;
		      if (keepOriginal)
			{
			  LanczosVectorStorage[MinPos].Copy(vec);
			}
		      else
			{
			  RealVector TmpV(vec);
			  vec=LanczosVectorStorage[MinPos];
			  LanczosVectorStorage[MinPos]=TmpV;
			}
		      VectorFlag |= (MinPos & VectorIndexMask);
		      VectorFlag |= SavedInMemory;
		      ProjectorLanczosVectorFlags[index]=VectorFlag;
		    }
		  else // no internal vectors stored in RAM -> write to disk, then.
		    {
		      sprintf(TmpOutputName,"vector-int.%d",index);
		      //cout << "Writing to disk: "<<TmpOutputName<<endl;
		      vec.WriteVector(TmpOutputName);
		      VectorFlag |= SavedOnDisk;
		      ProjectorLanczosVectorFlags[index] = VectorFlag;
		    }
		}
	    }
	}
      else // no memory storage at all, write things to disk straight away
	{
	  sprintf(TmpOutputName,"vector-int.%d",index);
	  vec.WriteVector(TmpOutputName);
	  //cout << "Writing to disk: "<<TmpOutputName<<endl;
	  VectorFlag |= SavedOnDisk;
	  ProjectorLanczosVectorFlags[index] = VectorFlag;
	}
    }
  return true;
}

// reread vector
// vec = vector to be retrieved
// index = vector index in Lanczos routine
// mainLanczos = flag indicating whether vector is part of main lanczos algorithm
// keepCopy = flag indicating whether the saved vector still needs to be kept in memory after reloading
void ProjectedReorthogonalizedLanczosAlgorithmDiskStorage::ReadVector(RealVector &vec, int index, bool mainLanczos, bool keepCopy)
{
  int VectorFlags;
  if (mainLanczos == true)
    VectorFlags = MainLanczosVectorFlags[index];
  else
    VectorFlags = ProjectorLanczosVectorFlags[index];
  if (VectorFlags==0)
    {
      cout << "Error: cannot retrieve vector "<<index<<" for "<<(mainLanczos?"main":"internal")<<" Lanczos!"<<endl;
      exit(-1);
    }
  if (VectorFlags & SavedInMemory)
    {
      int StoragePos = (VectorFlags & VectorIndexMask);
      if (keepCopy)
	{
	  vec.Copy(LanczosVectorStorage[StoragePos]);
	}
      else
	{
	  RealVector TmpV(vec);
	  vec=LanczosVectorStorage[StoragePos];
	  LanczosVectorStorage[StoragePos]=TmpV;
	  if (LanczosVectorStorage[StoragePos].GetVectorDimension()==0) // in case vector was deleted, recreate space
	    LanczosVectorStorage[StoragePos].Resize(this->Hamiltonian->GetHilbertSpaceDimension());
	  // note changes in memory
	  VectorFlags &= ~SavedInMemory;
	  VectorFlags &= ~VectorIndexMask;
	  if (mainLanczos == true)
	    MainLanczosVectorFlags[index] = VectorFlags;
	  else
	    ProjectorLanczosVectorFlags[index] = VectorFlags;
	}
    }
  else if (VectorFlags & SavedOnDisk)
    {
      if (mainLanczos == true)
	sprintf(TmpOutputName, "vector.%d", index);
      else
	sprintf(TmpOutputName, "vector-int.%d", index);
      vec.ReadVector(TmpOutputName);
    }
  else if (VectorFlags & SavedOnSwap1)
    {
      sprintf(TmpOutputName, "vector.swap-1");
      vec.ReadVector(TmpOutputName);
    }
  else if (VectorFlags & SavedOnSwap2)
    {
      sprintf(TmpOutputName, "vector.swap-2");
      vec.ReadVector(TmpOutputName);
    }
  return;
}

