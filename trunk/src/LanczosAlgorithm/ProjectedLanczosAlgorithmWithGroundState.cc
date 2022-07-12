///////////////////////////////////////////////////////////////////////////////
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


#include "LanczosAlgorithm/ProjectedLanczosAlgorithmWithGroundState.h"
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
#include <sys/time.h>
#include <iostream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;



// default constructor
//
// projectors = operators to use for projection after each application of the Hamiltonian
// nbrProjectors = number of separate projectors
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration
// nbrStorageVectors = number of vectors that can be held in memory (in addition to minimum of 3)
// diskFlag = use disk storage to increase speed of ground state calculation
// resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
// projectorPrecision = precision required for projector operations
// restartProjection = flag indicating whether projection should be restarted in precision not reached
ProjectedLanczosAlgorithmWithGroundState::ProjectedLanczosAlgorithmWithGroundState(AbstractHamiltonian** projectors, int nbrProjectors, AbstractArchitecture* architecture, int maxIter, int nbrStorageVectors, int projectorIterMax, bool diskFlag, bool resumeDiskFlag, double projectorPrecision, bool restartProjection)
{
  this->Flag.Initialize();
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
      MainIterMax=maxIter;
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(maxIter, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(maxIter);
    }
  else
    {
      MainIterMax=1<<12;
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
    }
  this->Architecture = architecture;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->NbrEigenvalue = 1;

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
  TmpOutputName = new char[50];
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ProjectedLanczosAlgorithmWithGroundState::ProjectedLanczosAlgorithmWithGroundState(const ProjectedLanczosAlgorithmWithGroundState& algorithm) 
{
  this->Flag = algorithm.Flag;
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

  TmpOutputName = new char[50];
}

// destructor
//

ProjectedLanczosAlgorithmWithGroundState::~ProjectedLanczosAlgorithmWithGroundState() 
{
  if (this->OrthogonalizationSet != 0)
    delete[] this->OrthogonalizationSet;
  if (this->OrthogonalizationSetFileNames != 0)
    {
      for (int i = 0; i < this->OrthogonalizationSetSize; ++i)
	delete[] this->OrthogonalizationSetFileNames[i];
      delete[] this->OrthogonalizationSetFileNames;
    }
  // clean up core of Projector
  this->ClearProjectorLanczosAlgorithm();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] this->VectorStorageFlags;
      delete [] this->MainLanczosVectorFlags;
      delete [] this->ProjectorLanczosVectorFlags;
      delete [] this->LanczosVectorStorage;
    }
  delete [] TmpOutputName;
}

// initialize Lanczos algorithm with a random vector
//

void ProjectedLanczosAlgorithmWithGroundState::InitializeLanczosAlgorithm() 
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
      this->ExternalOrthogonalization(this->V1);
      this->V1 /= this->V1.Norm();
      // project initial vector
      for (int p=0; p<NbrProjectors; ++p)
	this->ProjectVector(p);
      if (this->DiskFlag == false)
	this->InitialState = RealVector (this->V1, true);
      else
	{
	  // cout << "calling this->SaveVector(V1,0) on line "<<__LINE__<<endl;
	  this->SaveVector(V1,0,true,true);
	}
      this->Index = 0;
      this->TridiagonalizedMatrix.Resize(0, 0);
    }
  else
    {
      this->ReadState();
    }
  for (int i=0; i<NbrStorageVectors; ++i)
    this->LanczosVectorStorage[i].Resize(Dimension);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ProjectedLanczosAlgorithmWithGroundState::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  if (this->ResumeDiskFlag == false)
    {
      this->V1 = vector;
      this->V2 = RealVector (Dimension);
      this->V3 = RealVector (Dimension);
      if (this->OrthogonalizationSetSize > 0)
	{
	  this->ExternalOrthogonalization(this->V1);
	  this->V1 /= this->V1.Norm();
	}
      // project initial vector
      bool RequireReorthogonalization = false;
      for (int p=0; p<NbrProjectors; ++p)
	RequireReorthogonalization |= this->ProjectVector(p);
      if (this->DiskFlag == false)
	this->InitialState = RealVector (vector, true);
      else
	{
	  //cout << "calling this->SaveVector(V1,0,true,true) on line "<<__LINE__<<endl;
	  this->SaveVector(V1,0,true,true);
	}
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
  for (int i=0; i<NbrStorageVectors; ++i)
    this->LanczosVectorStorage[i].Resize(Dimension);
}

// force orthogonalization with respect to a set of vectors
//
// fileName = name of the file describing the set of vectors
// return value = true if no error occured

bool ProjectedLanczosAlgorithmWithGroundState::ForceOrthogonalization(char* fileName)
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

Vector* ProjectedLanczosAlgorithmWithGroundState::GetEigenstates(int nbrEigenstates)
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

Vector& ProjectedLanczosAlgorithmWithGroundState::GetGroundState()
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
	  this->ExternalOrthogonalization(this->V3);
	  this->V3 /= this->V3.Norm();

	  // apply projector here
	  // (V2 already saved as initial vector)
	  // swap V1, V3
	  {
	    RealVector TmpV (this->V1);
	    this->V1 = this->V3;
	    this->V3 = TmpV; 
	  }
	  //cout << "Line "<<__LINE__<<": V1="<<endl<<V1;
	  bool RequireReorthogonalization = false;
	  for (int p=0; p<NbrProjectors; ++p)
	    RequireReorthogonalization |= this->ProjectVector(p);
	  // swap V1, V2
	  {
	    RealVector TmpV (this->V1);
	    this->V1 = this->V3;
	    this->V3 = TmpV; 
	  }
	  this->V2.Copy(this->InitialState);
	  
	  // orthogonalizing again, if needed
	  if (RequireReorthogonalization)
	    {
	      double NewScalarProd = this->V2 * this->V3;
	      this->V3.AddLinearCombination(-NewScalarProd, this->V2);
	      this->ExternalOrthogonalization(this->V3);
	      this->V3 /= this->V3.Norm();
	    }
	  this->Swap1Index = -1;
	  this->Swap2Index = -1;
	  this->SaveVector(V3,1,true,true);		  
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
	      this->ExternalOrthogonalization(this->V1);
	      this->V1 /= this->V1.Norm();

	      bool RequireReorthogonalization = false;
	      for (int p=0; p<NbrProjectors; ++p)
		RequireReorthogonalization |= this->ProjectVector(p);
	      this->ReadVector(V2,i-2);
	      this->ReadVector(V3,i-1);
	      // end projector
	  
	      // reorthogonalize once more, if needed:
	      if (RequireReorthogonalization)
		{
		  // recalculate scalar products
		  RealVector* TmpVectorScalarProduct[2];
		  double TmpScalarProduct[2];
		  TmpVectorScalarProduct[0] = &(this->V1);
		  TmpVectorScalarProduct[1] = &(this->V2);
		  MultipleRealScalarProductOperation Operation (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
		  Operation.ApplyOperation(this->Architecture);
		  // perform subtractions
		  TmpVector[0] = this->V1;
		  TmpVector[1] = this->V2;
		  TmpCoefficient[0] = -TmpScalarProduct[0];
		  TmpCoefficient[1] = -TmpScalarProduct[1];
		  AddRealLinearCombinationOperation Operation1 (&(this->V3),  TmpVector, 2, TmpCoefficient);
		  Operation1.ApplyOperation(this->Architecture);	  
		  this->ExternalOrthogonalization(this->V3);
		  this->V3 /= this->V3.Norm();
		}
	      delete[] TmpVector;
	  
	      this->SaveVector(V1, i, true, true);
	      
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
	  this->ReadVector(V1,0);
	  this->GroundState.Copy(this->V1, TmpComponents[0]);
	  for (int i = 1; i < this->DiagonalizedMatrix.GetNbrRow(); ++i)
	    {
	      this->ReadVector(V1,i);
	      this->GroundState.AddLinearCombination(TmpComponents[i], this->V1);	      
	      cout << i << "/" << this->DiagonalizedMatrix.GetNbrRow() << "           \r";
	      cout.flush();
	    }	  
	}
      cout << endl;
      this->ExternalOrthogonalization(this->GroundState);
      this->GroundState /= this->GroundState.Norm();
      this->GroundStateFlag = true;
      delete[] TmpComponents;
    }
  return this->GroundState;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void ProjectedLanczosAlgorithmWithGroundState::RunLanczosAlgorithm (int nbrIter) 
{
  this->GroundStateFlag = false;
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
      this->ExternalOrthogonalization(this->V2);
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
      this->ReadVector(V1,0);      
      
      // end projector

      // orthogonalizing again, if needed
      if (RequireReorthogonalization)
	{
	  double NewScalarProd = this->V1 * this->V2;
	  this->V2.AddLinearCombination(-NewScalarProd, this->V1);
	  this->ExternalOrthogonalization(this->V2);
	  this->V2 /= this->V2.Norm();
	}
      
      if (this->DiskFlag == true)
	{
	  // cout << "calling this->SaveVector(V2,1) on line "<<__LINE__<<endl;
	  this->SaveVector(V2,1,true,true);
	}
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
	  this->ExternalOrthogonalization(this->V3);
	  this->V3 /= this->V3.Norm();

	  // apply projector here
	  // vector0..i-1 are already saved
	  if (this->DiskFlag == false) // have already saved contents of V1 above if DiskFlag is true
	    {
	      // cout << "calling this->SaveVector(V2,"<<i<<") on line "<<__LINE__<<endl;
	      this->SaveVector(V2, i-1); // do not need copy here, once more
	    }

	  // swap V1, V3
	  {
	    RealVector TmpV (this->V1);
	    this->V1 = this->V3;
	    this->V3 = TmpV;
	  }
	  bool RequireReorthogonalization = false;
	  for (int p=0; p<NbrProjectors; ++p)
	    RequireReorthogonalization |= this->ProjectVector(p);
	  // swap V1, V3
	  {
	    RealVector TmpV (this->V1);
	    this->V1 = this->V3;
	    this->V3 = TmpV; 
	  }
	  this->ReadVector(V1,i-2);
	  this->ReadVector(V2,i-1);
	  // end projector
	  
	  // reorthogonalize once more, if needed:
	  if (RequireReorthogonalization)
	    {
	      // recalculate scalar products
	      RealVector* TmpVectorScalarProduct[2];
	      TmpVectorScalarProduct[0] = &(this->V1);
	      TmpVectorScalarProduct[1] = &(this->V2);
	      MultipleRealScalarProductOperation Operation (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
	      Operation.ApplyOperation(this->Architecture);
	      // perform subtractions
	      TmpVector[0] = this->V1;
	      TmpVector[1] = this->V2;
	      TmpCoefficient[0] = -TmpScalarProduct[0];
	      TmpCoefficient[1] = -TmpScalarProduct[1];
	      AddRealLinearCombinationOperation Operation1 (&(this->V3),  TmpVector, 2, TmpCoefficient);
	      Operation1.ApplyOperation(this->Architecture);	  
	      this->ExternalOrthogonalization(this->V3);
	      this->V3 /= this->V3.Norm();
	    }
	  delete[] TmpVector;
	  
	  if (this->DiskFlag == true)
	    {
	      // cout << "calling this->SaveVector(V3,"<<i<<") on line "<<__LINE__<<endl;
	      this->SaveVector(V3, i, true, true);
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
	  this->ReadVector(V1,i-1);
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

bool ProjectedLanczosAlgorithmWithGroundState::TestConvergence ()
{
  if ((fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
       (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1)))))
    return true;
  else
    return false;
}


// project the vector stored in V1 to the groundstate of the given Projector number
// nbrProjector = number of projector to use for this run
// return = true if Lanczos was run, or false if no projection needed to be done
bool ProjectedLanczosAlgorithmWithGroundState::ProjectVector(int nbrProjector)
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

bool ProjectedLanczosAlgorithmWithGroundState::TestProjectorConvergence ()
{
  if ( (fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0)) < this->ProjectorPrecision) ||
       (fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0) - this->PreviousProjectorGroundstate) < 
	(this->ProjectorPrecision * fabs(this->InternalDiagonalizedMatrix.DiagonalElement(0)))))
    return true;
  else
    return false;
}


// run current internal Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrProjector = number of projector to use for this run
// nbrIter = number of iteration to do 

bool ProjectedLanczosAlgorithmWithGroundState::RunProjectorLanczosAlgorithm (int nbrProjector, int nbrIter) 
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
      this->SaveVector(V1, 0, false, true); // need to keep V1 in process
      this->V2.AddLinearCombination(-this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex), 
				    this->V1);
      this->ExternalOrthogonalization(this->V2);
      this->V2 /= this->V2.Norm();
      // cout << "calling this->SaveVector(V2,"<<1<<",false,true) on line "<<__LINE__<<endl;
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
  double* TmpCoefficient = new double[2];
  double TmpScalarProduct[2];
  for (int i = this->InternalIndex + 2; i < Dimension; i++)
    {
      RealVector* TmpVector = new RealVector[2];

      TmpVector[0] = this->V1;
      TmpVector[1] = this->V2;
      TmpCoefficient[0] = -this->InternalTridiagonalizedMatrix.UpperDiagonalElement(this->InternalIndex);
      TmpCoefficient[1] = -this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex + 1);
      AddRealLinearCombinationOperation Operation4 (&(this->V3),  TmpVector, 2, TmpCoefficient);
      Operation4.ApplyOperation(this->Architecture);
      delete[] TmpVector;
      this->ExternalOrthogonalization(this->V3);
      this->V3 /= this->V3.Norm();
      // cout << "calling this->SaveVector(V3,"<<i<<",false) on line "<<__LINE__<<endl;
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
      RealVector* TmpVectorScalarProduct[2];
      TmpVectorScalarProduct[0] = &(this->V1);
      TmpVectorScalarProduct[1] = &(this->V2);
      MultipleRealScalarProductOperation Operation2 (&(this->V3), TmpVectorScalarProduct, 2, TmpScalarProduct);
      Operation2.ApplyOperation(this->Architecture);
      this->InternalTridiagonalizedMatrix.UpperDiagonalElement(this->InternalIndex) = TmpScalarProduct[0];
      this->InternalTridiagonalizedMatrix.DiagonalElement(this->InternalIndex + 1) = TmpScalarProduct[1];
    }
  delete[] TmpCoefficient;
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

Vector& ProjectedLanczosAlgorithmWithGroundState::GetProjectorGroundState()
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
  this->ExternalOrthogonalization(this->V1);
  this->V1 /= this->V1.Norm();  
  delete[] TmpComponents;
  return this->V1;
}


// write current Lanczos state on disk
//
// return value = true if no error occurs

bool ProjectedLanczosAlgorithmWithGroundState::WriteState()
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

bool ProjectedLanczosAlgorithmWithGroundState::ReadState()
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
  sprintf(TmpOutputName, "vector.%d", this->Index);
  this->V1.ReadVector(TmpOutputName);
  sprintf(TmpOutputName, "vector.%d", (this->Index + 1));
  this->V2.ReadVector(TmpOutputName);
  sprintf(TmpOutputName, "vector.%d", (this->Index + 2));
  this->V3.ReadVector(TmpOutputName);
  return true;
}


// orthogonalize a vector with respect to a set of external vectors
//
// inputVector = reference on the vector whose component on the external set has to be removed 

void  ProjectedLanczosAlgorithmWithGroundState::ExternalOrthogonalization(RealVector& inputVector)
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


// initialize Lanczos algorithm starting from the vector currently held in V1
//
void ProjectedLanczosAlgorithmWithGroundState::InitializeProjectorLanczosAlgorithm()
{
  if (this->InternalIndex>0)
    this->ClearProjectorLanczosAlgorithm();
  // normalization and external orthogonalization should already have been done before calling this...
  // this->ExternalOrthogonalization(this->V1);
  // this->V1 /= this->V1.Norm();  
  this->InternalTridiagonalizedMatrix.ResizeAndClean(0, 0);
}


// clear storage associated with internal Lanczos algorithm
//
void ProjectedLanczosAlgorithmWithGroundState::ClearProjectorLanczosAlgorithm()
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
void ProjectedLanczosAlgorithmWithGroundState::ProjectorDiagonalize ()
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
bool ProjectedLanczosAlgorithmWithGroundState::SaveVectorTest(RealVector &vec, int index, bool mainLanczos, bool keepOriginal)
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
void ProjectedLanczosAlgorithmWithGroundState::ReadVectorTest(RealVector &vec, int index, bool mainLanczos, bool keepCopy)
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
bool ProjectedLanczosAlgorithmWithGroundState::SaveVector(RealVector &vec, int index, bool mainLanczos, bool keepOriginal)
{
  int StorageFlag=index & StorageIndexMask;
  int VectorFlag=0;
  
  if (mainLanczos==true)
    {
      StorageFlag |= MainLanczosVector;
      if (this->DiskFlag==true)
	{	  
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
			      int TmpIndex = (VectorStorageFlags[i]&VectorIndexMask);
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
			  cout << "Problem in ProjectedLanczosAlgorithmWithGroundState::SaveVector - should have saved Main vector no. MinIndex, but not flagged"<<endl;
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
      else
	{
	  // saving only up to two vectors of main Lanczos that are swapped out for internal lanczos run
	  if (NbrStorageVectors>2) // can save them in memory?
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
			  }
			else
			  {
			    RealVector TmpV(vec);
			    vec=LanczosVectorStorage[i];
			    LanczosVectorStorage[i]=TmpV;
			  }
			cout << "attention, overwriting vector "<<index<<" in RAM"<<endl;
			break;
		      }
		}
	      else
		{
		  // search for additional storage in live memory
		  if (NbrStoredMain<2) // have more space to give away?
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
			      int TmpIndex = (VectorStorageFlags[i]&VectorIndexMask);
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
		  else // have no more space to give away: overwrite oldest vector of main lanczos
		    {
		      int MinPos=0;
		      int MinIndex=MainIterMax;
		      for (int i=0; i<NbrStorageVectors; ++i)
			{
			  if ((VectorStorageFlags[i]&MainLanczosVector) &&
			      ((VectorStorageFlags[i]&StorageIndexMask) < MinIndex))
			    {
			      MinPos=i;
			      MinIndex=(VectorStorageFlags[i]&StorageIndexMask);
			    }
			}
		      // mark state as forgotten from memory
		      MainLanczosVectorFlags[index]&= ~SavedInMemory;
		      MainLanczosVectorFlags[index]&= ~VectorIndexMask;
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
	  else // save the state on disk, using only two vector files vector.swap-1 / vector.swap-2
	    {
	      bool UseSwap2;
	      if (index == Swap2Index)
		UseSwap2 = true;
	      else if (index == Swap1Index)
		UseSwap2 = false;
	      else UseSwap2 = (Swap2Index<Swap1Index);
	      if (UseSwap2)
		{
		  sprintf(TmpOutputName,"vector.swap-2");
		  //cout << "Writing to disk: "<<TmpOutputName<<endl;
		  vec.WriteVector(TmpOutputName);
		  VectorFlag |= SavedOnSwap2;
		  MainLanczosVectorFlags[index] = VectorFlag;
		  if (Swap2Index>-1)
		    MainLanczosVectorFlags[Swap2Index] = 0;
		  Swap2Index = index;
		}
	      else
		{
		  sprintf(TmpOutputName,"vector.swap-1");
		  //cout << "Writing to disk: "<<TmpOutputName<<endl;
		  vec.WriteVector(TmpOutputName);
		  VectorFlag |= SavedOnSwap1;
		  MainLanczosVectorFlags[index] = VectorFlag;
		  if (Swap1Index>-1)
		    MainLanczosVectorFlags[Swap1Index] = 0;
		  Swap1Index = index;
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
		      cout << "Problem in ProjectedLanczosAlgorithmWithGroundState::SaveVector - should have saved Main vector no. MinIndex, but not flagged"<<endl;
		      exit(-1);
		    }
		  // mark state as forgotten
		  MainLanczosVectorFlags[MinIndex]&= ~SavedInMemory;
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
void ProjectedLanczosAlgorithmWithGroundState::ReadVector(RealVector &vec, int index, bool mainLanczos, bool keepCopy)
{
  int VectorFlags;
  if (mainLanczos == true)
    VectorFlags = MainLanczosVectorFlags[index];
  else
    VectorFlags = ProjectorLanczosVectorFlags[index];
  if (VectorFlags==0)
    {
      if ((index==0) && (this->DiskFlag==false))
	{
	  vec.Copy(this->InitialState);
	  return;
	}
      else
	{
	  cout << "Error: cannot retrieve vector "<<index<<" for "<<(mainLanczos?"main":"internal")<<" Lanczos!"<<endl;
	  exit(-1);
	}
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
