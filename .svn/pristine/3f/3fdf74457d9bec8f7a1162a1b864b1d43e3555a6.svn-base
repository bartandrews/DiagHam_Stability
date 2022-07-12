////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of operations that compute the matrix elements            //
//                  of the many-body interaction on the torus                 // 
//                                                                            //
//                        last modification : 06/03/2015                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusComputeMatrixElementOperation.h"

#include <algorithm>


// constructor 
//
// hamiltonian = pointer to the generic n-body Hamiltonian
// nbrUniqueMatrixElements = number of unique matrix elements
// momentumSectorIndices = array that contains the momentum sector of each matrix element
// j1Indices = array that contains the creation indices of each matrix element
// j2Indices = array that contains the annihilation indices of each matrix element
// matrixElements = array where the matrix elements will be stored
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHETorusComputeMatrixElementOperation::FQHETorusComputeMatrixElementOperation(ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian* hamiltonian, long nbrUniqueMatrixElements,
									       int* momentumSectorIndices, int* j1Indices, int* j2Indices, Complex* matrixElements,
									       int nbrMPIStage, int nbrSMPStage)
{
  this->Hamiltonian = hamiltonian;
  this->TwistedHamiltonian = 0;
  this->NbrUniqueMatrixElements = nbrUniqueMatrixElements;
  this->MomentumSectorIndices = momentumSectorIndices;
  this->J1Indices = j1Indices;
  this->J2Indices = j2Indices;
  this->MatrixElements = matrixElements;

  this->FirstComponent = 0;
  this->NbrComponent = (int) nbrUniqueMatrixElements;
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = nbrUniqueMatrixElements;
  this->OperationType = AbstractArchitectureOperation::FQHETorusComputeMatrixElementOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}
    
// constructor 
//
// hamiltonian = pointer to the generic n-body Hamiltonian
// nbrUniqueMatrixElements = number of unique matrix elements
// momentumSectorIndices = array that contains the momentum sector of each matrix element
// j1Indices = array that contains the creation indices of each matrix element
// j2Indices = array that contains the annihilation indices of each matrix element
// matrixElements = array where the matrix elements will be stored
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHETorusComputeMatrixElementOperation::FQHETorusComputeMatrixElementOperation(ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian* hamiltonian, long nbrUniqueMatrixElements,
									       int* momentumSectorIndices, int* j1Indices, int* j2Indices, Complex* matrixElements,
									       int nbrMPIStage, int nbrSMPStage)
{
  this->Hamiltonian = 0;
  this->TwistedHamiltonian = hamiltonian;
  this->NbrUniqueMatrixElements = nbrUniqueMatrixElements;
  this->MomentumSectorIndices = momentumSectorIndices;
  this->J1Indices = j1Indices;
  this->J2Indices = j2Indices;
  this->MatrixElements = matrixElements;

  this->FirstComponent = 0;
  this->NbrComponent = (int) nbrUniqueMatrixElements;
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = nbrUniqueMatrixElements;
  this->OperationType = AbstractArchitectureOperation::FQHETorusComputeMatrixElementOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}
    
// copy constructor 
//
// operation = reference on operation to copy

FQHETorusComputeMatrixElementOperation::FQHETorusComputeMatrixElementOperation(const FQHETorusComputeMatrixElementOperation & operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHETorusComputeMatrixElementOperation;
  this->Hamiltonian = operation.Hamiltonian;
  this->TwistedHamiltonian = operation.TwistedHamiltonian;
  this->NbrUniqueMatrixElements = operation.NbrUniqueMatrixElements;
  this->MomentumSectorIndices = operation.MomentumSectorIndices;
  this->J1Indices = operation.J1Indices;
  this->J2Indices = operation.J2Indices;
  this->MatrixElements = operation.MatrixElements;
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;
}
  
// destructor
//

FQHETorusComputeMatrixElementOperation::~FQHETorusComputeMatrixElementOperation()
{
}
  
// set range of indices
//
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHETorusComputeMatrixElementOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->LargeFirstComponent = firstComponent;
  this->LargeNbrComponent = nbrComponent;
  this->FirstComponent = (int) firstComponent;
  this->NbrComponent = (int) nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHETorusComputeMatrixElementOperation::Clone()
{
  return new FQHETorusComputeMatrixElementOperation (*this);
}

// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::RawApplyOperation()
{
  int LastIndex = this->FirstComponent + this->NbrComponent;
  int TmpNBodyValue = 0;
  if (this->Hamiltonian != 0)
    {
      TmpNBodyValue = this->Hamiltonian->NBodyValue;
    }
  else
    {
      TmpNBodyValue = this->TwistedHamiltonian->NBodyValue;
    }
int* TmpNIndices =  new int [TmpNBodyValue];
  int* TmpMIndices =  new int [TmpNBodyValue];
  double* QxValues = new double [TmpNBodyValue];
  double* QyValues = new double [TmpNBodyValue];
  double* Q2Values = new double [TmpNBodyValue];
  double* CosineCoefficients = new double [TmpNBodyValue];
  if (((this->Hamiltonian != 0) && (this->Hamiltonian->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic))
      || ((this->TwistedHamiltonian != 0) && (this->TwistedHamiltonian->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)))
    {
      int NbrPermutations = 1;
      for (int i = 2; i <= TmpNBodyValue; ++i)
	NbrPermutations *= i;
      int** Permutations = new int*[NbrPermutations]; 
      double* PermutationSign = new double[NbrPermutations]; 
      Permutations[0] = new int [TmpNBodyValue];
      for (int i = 0; i < TmpNBodyValue; ++i)
	Permutations[0][i] = i;
      PermutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  Permutations[i] = new int [TmpNBodyValue];
	  for (int j = 0; j < TmpNBodyValue; ++j)
	    Permutations[i][j] = Permutations[i - 1][j];
	  int* TmpArrayPerm = Permutations[i];
	  int Pos1 = TmpNBodyValue - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = TmpNBodyValue - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = TmpNBodyValue - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  PermutationSign[i] = TmpSign;
	}
      
      if (this->Hamiltonian != 0)
	{
	  for (int Index = this->FirstComponent;  Index < LastIndex; ++Index)
	    {
	      int J1 = this->J1Indices[Index];
	      int J2 = this->J2Indices[Index];
	      int MomentumSector = this->MomentumSectorIndices[Index];
	      double TmpInteraction = 0.0;
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm1 = Permutations[l1];
		  for (int k = 0; k < TmpNBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->Hamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J1 * TmpNBodyValue) + TmpPerm1[k]];
		    }
		  for (int l2 = 0; l2 < NbrPermutations; ++l2)
		    {
		      int* TmpPerm2 = Permutations[l2];
		      for (int k = 0; k < TmpNBodyValue; ++k)
			{
			  TmpMIndices[k] = this->Hamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J2 * TmpNBodyValue) + TmpPerm2[k]];
			}
		      double Tmp = this->Hamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
										     QxValues, QyValues, Q2Values, CosineCoefficients);
		      TmpInteraction += PermutationSign[l1] * PermutationSign[l2] * Tmp;
		    }
		}
	      this->MatrixElements[Index] = TmpInteraction;
	    }
	}
      else
	{
	  for (int Index = this->FirstComponent;  Index < LastIndex; ++Index)
	    {
	      int J1 = this->J1Indices[Index];
	      int J2 = this->J2Indices[Index];
	      int MomentumSector = this->MomentumSectorIndices[Index];
	      Complex TmpInteraction = 0.0;
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm1 = Permutations[l1];
		  for (int k = 0; k < TmpNBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->TwistedHamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J1 * TmpNBodyValue) + TmpPerm1[k]];
		    }
		  for (int l2 = 0; l2 < NbrPermutations; ++l2)
		    {
		      int* TmpPerm2 = Permutations[l2];
		      for (int k = 0; k < TmpNBodyValue; ++k)
			{
			  TmpMIndices[k] = this->TwistedHamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J2 * TmpNBodyValue) + TmpPerm2[k]];
			}
		      Complex Tmp = this->TwistedHamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
											     QxValues, QyValues, Q2Values, CosineCoefficients);
		      TmpInteraction += PermutationSign[l1] * PermutationSign[l2] * Tmp;
		    }
		}
	      this->MatrixElements[Index] = TmpInteraction;
	    }
	}
      delete[] PermutationSign;
      for (int i = 0; i < NbrPermutations; ++i)
	{
	  delete[] Permutations[i];
	}
      delete[] Permutations;
    }
  else
    {
      if (this->Hamiltonian != 0)
	{
	}
      else
	{
	  for (int Index = this->FirstComponent;  Index < LastIndex; ++Index)
	    {
	      int J1 = this->J1Indices[Index];
	      int J2 = this->J2Indices[Index];
	      int MomentumSector = this->MomentumSectorIndices[Index];
	      for (int k = 0; k < TmpNBodyValue; ++k)
		{
		  TmpNIndices[k] = this->TwistedHamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J1 * TmpNBodyValue) + k];
		  TmpMIndices[k] = this->TwistedHamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J2 * TmpNBodyValue) + k];
		}
	      Complex TmpInteraction = this->TwistedHamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
												QxValues, QyValues, Q2Values, CosineCoefficients);
	      while (std::prev_permutation(TmpMIndices, TmpMIndices  + TmpNBodyValue))
		{
		  TmpInteraction += this->TwistedHamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, QxValues, QyValues, Q2Values, CosineCoefficients);		    
		}
	      while (std::prev_permutation(TmpNIndices, TmpNIndices  + TmpNBodyValue))
		{
		  for (int k = 0; k < TmpNBodyValue; ++k)
		    {
		      TmpMIndices[k] = this->TwistedHamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J2 * TmpNBodyValue) + k];
		    }
		  TmpInteraction += this->TwistedHamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, QxValues, QyValues, Q2Values, CosineCoefficients);
		  while (std::prev_permutation(TmpMIndices, TmpMIndices  + TmpNBodyValue))
		    {
		      TmpInteraction += this->TwistedHamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, QxValues, QyValues, Q2Values, CosineCoefficients);
		    }
		}
	      this->MatrixElements[Index] = TmpInteraction;
	    }
	}
    }
  delete[] TmpNIndices;
  delete[] TmpMIndices;
  delete[] QxValues;
  delete[] QyValues;
  delete[] Q2Values;
  delete[] CosineCoefficients;
  return true;
}

// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{
  int TmpNbrComponents = this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  
  int NbrStages = this->NbrSMPStage * architecture->GetNbrThreads();
  int StageIndex = 0;
   bool LockFlag = false;
   while (StageIndex < NbrStages) 
     {
       if (LockFlag == false)   
         {
           architecture->LockMutex();
           LockFlag = true;
         }
       StageIndex = this->SMPStages[0];
       if (StageIndex < NbrStages) 
         {         
           this->SMPStages[0]++;   
           architecture->UnLockMutex();
           LockFlag = false;
           this->SetIndicesRange(TmpFirstComponent + this->GetRankChunkStart(TmpNbrComponents, StageIndex,  NbrStages),  
				 this->GetRankChunkSize(TmpNbrComponents, StageIndex,  NbrStages));
           this->RawApplyOperation();
           ++StageIndex;
         }
     }
   if (LockFlag == true)
     {
       architecture->UnLockMutex();
       LockFlag = false;
     }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  FQHETorusComputeMatrixElementOperation** TmpOperations = new FQHETorusComputeMatrixElementOperation* [architecture->GetNbrThreads()];
  for(int i = 0; i < architecture->GetNbrThreads(); ++i)
    TmpOperations[i] = (FQHETorusComputeMatrixElementOperation *) this->Clone();
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  return true;
//   int Step = this->NbrComponent / architecture->GetNbrThreads();
//   if (Step == 0)
//     Step = this->NbrComponent;
//   this->SMPStages[0] = 0;
//   int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
//   int TmpFirstComponent = this->FirstComponent;
//   int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
//   FQHETorusComputeMatrixElementOperation** TmpOperations = new  FQHETorusComputeMatrixElementOperation* [architecture->GetNbrThreads()];
  
   
//   for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
//     {
//       TmpOperations[i] = (FQHETorusComputeMatrixElementOperation *) this->Clone();
//       architecture->SetThreadOperation(TmpOperations[i], i);
//     }
  
//   architecture->SendJobsRoundRobin();

//   for (int i = 0; i < architecture->GetNbrThreads(); i++)
//     {
//       delete TmpOperations[i];
//     }
//   delete[] TmpOperations;
//   return true;
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
//   int Step = this->NbrComponent / architecture->GetNbrNodes();
//   int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
//   int TmpNbrComponent = Step;
//   if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
//     {
//       TmpNbrComponent += this->NbrComponent % Step;
//     }
//   this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
//   switch (architecture->GetArchitectureID())
//     {	 
//     case AbstractArchitecture::MixedMPISMP:
//       this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
//       break;
//     default:
//       this->RawApplyOperation();
//       break;
//     }		
//   MPI::COMM_WORLD.Barrier();
//   if (this->OutputState != 0)
//     architecture->SumVector(*(this->OutputState));	
//   else
//     architecture->SumVector(*(this->ComplexOutputState));	
  return true;
#else
  return this->RawApplyOperation();
#endif
}
  
