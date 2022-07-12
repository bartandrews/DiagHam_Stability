////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of spin chain multiple entanglement spectrum calculation        //
//                                                                            //
//                        last modification : 06/08/2016                      //
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
#include "Architecture/ArchitectureOperation/SpinChainMultipleEntanglementSpectrumOperation.h"
#include "Operator/AbstractOperator.h"
#include "GeneralTools/ArrayTools.h"


// constructor 
//
// space = pointer to the Hilbert space
// eigenstates = matrix that contains the eigenstates
// firstEigenstate = first eigenstate to consider
// lastEigenstate = last eigenstate to consider
// subsystemSize = number of sites in the subsystem
// subsystemSz = twice the total Sz value for the subsystem
// subsystemSites = array that describes the subsystem (if null, consider the subsystemSize first sites)

SpinChainMultipleEntanglementSpectrumOperation::SpinChainMultipleEntanglementSpectrumOperation(AbstractSpinChain* space, RealMatrix& eigenstates, int firstEigenstate, 
											       int lastEigenstate, int subsystemSize, int subsystemSz, int* subsystemSites)
{
  this->SubsystemSize = subsystemSize;
  this->SubsystemSz = subsystemSz;
  this->RealEigenstates = eigenstates;
  this->OperationType = AbstractArchitectureOperation::SpinChainMultipleEntanglementSpectrum;  
  this->FirstState = firstEigenstate;
  this->NbrStates = lastEigenstate - firstEigenstate + 1;
  this->EntanglementSpectrumDimension = 0;
  this->EntanglementSpectra = 0;
  this->Space = (AbstractSpinChain*) space->Clone();
  if (subsystemSites != 0)
    {
      this->SubsystemSites = new int[this->SubsystemSize];
      for (int i = 0; i < this->SubsystemSize; ++i)
	{
	  this->SubsystemSites[i] = subsystemSites[i];
	}
    }
  else
    {
      this->SubsystemSites = 0;
    }
}

// constructor 
//
// space = pointer to the Hilbert space
// eigenstates = matrix that contains the eigenstates
// firstEigenstate = first eigenstate to consider
// lastEigenstate = last eigenstate to consider
// subsystemSize = number of sites in the subsystem
// subsystemSz = twice the total Sz value for the subsystem
// subsystemSites = array that describes the subsystem (if null, consider the subsystemSize first sites)

SpinChainMultipleEntanglementSpectrumOperation::SpinChainMultipleEntanglementSpectrumOperation(AbstractSpinChain* space, ComplexMatrix& eigenstates, int firstEigenstate, 
											       int lastEigenstate, int subsystemSize, int subsystemSz, int* subsystemSites)
{
  this->SubsystemSize = subsystemSize;
  this->SubsystemSz = subsystemSz;
  this->ComplexEigenstates = eigenstates;
  this->OperationType = AbstractArchitectureOperation::SpinChainMultipleEntanglementSpectrum;  
  this->FirstState = firstEigenstate;
  this->NbrStates = lastEigenstate - firstEigenstate + 1;
  this->EntanglementSpectrumDimension = 0;
  this->EntanglementSpectra = 0;
  this->Space = (AbstractSpinChain*) space->Clone();
  if (subsystemSites != 0)
    {
      this->SubsystemSites = new int[this->SubsystemSize];
      for (int i = 0; i < this->SubsystemSize; ++i)
	{
	  this->SubsystemSites[i] = subsystemSites[i];
	}
    }
  else
    {
      this->SubsystemSites = 0;
    }
}

// copy constructor 
//
// operation = reference on operation to copy

SpinChainMultipleEntanglementSpectrumOperation::SpinChainMultipleEntanglementSpectrumOperation(const SpinChainMultipleEntanglementSpectrumOperation& operation)
{
  this->RealEigenstates = operation.RealEigenstates;
  this->ComplexEigenstates = operation.ComplexEigenstates;
  this->SubsystemSize = operation.SubsystemSize;
  this->SubsystemSz = operation.SubsystemSz;
  this->FirstState = operation.FirstState;
  this->NbrStates = operation.NbrStates;
  this->EntanglementSpectra = 0;
  this->EntanglementSpectrumDimension = operation.EntanglementSpectrumDimension;
  this->Space = (AbstractSpinChain*) (operation.Space->Clone());
  if (operation.SubsystemSites != 0)
    {
      this->SubsystemSites = new int[this->SubsystemSize];
      for (int i = 0; i < this->SubsystemSize; ++i)
	{
	  this->SubsystemSites[i] = operation.SubsystemSites[i];
	}
    }
  else
    {
      this->SubsystemSites = 0;
    }
}
  
// destructor
//

SpinChainMultipleEntanglementSpectrumOperation::~SpinChainMultipleEntanglementSpectrumOperation()
{
  if (this->SubsystemSites != 0)
    delete[] this->SubsystemSites;
  delete this->Space;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* SpinChainMultipleEntanglementSpectrumOperation::Clone()
{
  return new SpinChainMultipleEntanglementSpectrumOperation (*this);
}
  

// set the number of operators that have to be locally evaluated
// 
// firstState = index of the first state to evaluate 
// nbrStates =number of states that have to be evaluated bu the local architecture

void SpinChainMultipleEntanglementSpectrumOperation::SetStateRange(int firstState, int nbrStates)
{
  this->FirstState = firstState;
  this->NbrStates = nbrStates;
}

// apply operation(architecture independent)
//
// return value = true if no error occurs

bool SpinChainMultipleEntanglementSpectrumOperation::RawApplyOperation()
{
  int LastState = this->FirstState + this->NbrStates;
  this->EntanglementSpectra = new double*[this->NbrStates];
  if (this->RealEigenstates.GetNbrRow() != 0)
    {
      for (int Index = this->FirstState; Index < LastState; ++Index)
	{
	  RealMatrix PartialEntanglementMatrix;
 	  if (this->SubsystemSites == 0)
 	    PartialEntanglementMatrix = this->Space->EvaluatePartialEntanglementMatrix(this->SubsystemSize, this->SubsystemSz, this->RealEigenstates[Index]);
 	  else
 	    PartialEntanglementMatrix = this->Space->EvaluatePartialEntanglementMatrix(this->SubsystemSites, this->SubsystemSize, this->SubsystemSz, this->RealEigenstates[Index]);
	  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
	    {
	      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
	      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		}
	      for (int i = 0; i < TmpDimension; ++i)
		TmpValues[i] *= TmpValues[i];   
	      SortArrayDownOrdering(TmpValues, TmpDimension);
	      this->EntanglementSpectrumDimension = TmpDimension;
	      this->EntanglementSpectra[Index - this->FirstState] = TmpValues;
	    }
	  else
	    {
	      double* TmpValues = new double[1];
	      TmpValues[0] = 0.0;
	      if (PartialEntanglementMatrix.GetNbrRow() == 1)
		{
		  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
		    TmpValues[0] += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
		}
	      else
		{
		  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
		    TmpValues[0] += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
		}
	      this->EntanglementSpectrumDimension = 1;
	      this->EntanglementSpectra[Index - this->FirstState] = TmpValues;
	    }
	}
    }
  else
    {
      for (int Index = this->FirstState; Index < LastState; ++Index)
	{
	  ComplexMatrix PartialEntanglementMatrix;
	  if (this->SubsystemSites == 0)
	    PartialEntanglementMatrix = this->Space->EvaluatePartialEntanglementMatrix(this->SubsystemSize, this->SubsystemSz, this->ComplexEigenstates[Index]);
	  else
	    PartialEntanglementMatrix = this->Space->EvaluatePartialEntanglementMatrix(this->SubsystemSites, this->SubsystemSize, this->SubsystemSz, this->ComplexEigenstates[Index]);
	  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
	    {
	      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
	      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		}
	      for (int i = 0; i < TmpDimension; ++i)
		TmpValues[i] *= TmpValues[i];   
	      SortArrayDownOrdering(TmpValues, TmpDimension);
	      this->EntanglementSpectrumDimension = TmpDimension;
	      this->EntanglementSpectra[Index - this->FirstState] = TmpValues;
	    }
	  else
	    {
	      double* TmpValues = new double[1];
	      TmpValues[0] = 0.0;
	      if (PartialEntanglementMatrix.GetNbrRow() == 1)
		{
		  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
		    TmpValues[0] += SqrNorm(PartialEntanglementMatrix[i][0]);
		}
	      else
		{
	      for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
		TmpValues[0] += SqrNorm(PartialEntanglementMatrix[0][i]);				  
		}
	      this->EntanglementSpectrumDimension = 1;
	      this->EntanglementSpectra[Index - this->FirstState] = TmpValues;
	    }
	}
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool SpinChainMultipleEntanglementSpectrumOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrStates / architecture->GetNbrThreads();
  if (Step == 0)
    Step = 1;
  int TmpFirstState = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  SpinChainMultipleEntanglementSpectrumOperation** TmpOperations = new SpinChainMultipleEntanglementSpectrumOperation* [architecture->GetNbrThreads()];
  for (int OperatorIndex = 0; OperatorIndex < ReducedNbrThreads; ++OperatorIndex)
    {
      TmpOperations[OperatorIndex] = (SpinChainMultipleEntanglementSpectrumOperation*) this->Clone();
      if ((TmpFirstState + Step) >= this->NbrStates)
	{
	  TmpOperations[OperatorIndex]->SetStateRange(TmpFirstState, this->NbrStates - TmpFirstState);
	  Step = 0;
	  TmpFirstState = 0;
	}
      else
	{
	  TmpOperations[OperatorIndex]->SetStateRange(TmpFirstState,Step);
	  TmpFirstState += Step;
	}
      architecture->SetThreadOperation(TmpOperations[OperatorIndex], OperatorIndex);
    }
  TmpOperations[ReducedNbrThreads] = (SpinChainMultipleEntanglementSpectrumOperation*) this->Clone();
  if (Step > 0)
    TmpOperations[ReducedNbrThreads]->SetStateRange(TmpFirstState, this->NbrStates - TmpFirstState);
  else
    TmpOperations[ReducedNbrThreads]->SetStateRange(TmpFirstState, Step);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();

  this->EntanglementSpectra = new double*[this->NbrStates];
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      int TmpLastState = TmpOperations[i]->FirstState + TmpOperations[i]->NbrStates;
      for (int j = TmpOperations[i]->FirstState; j < TmpLastState; ++j)
	{
	  this->EntanglementSpectra[j - this->FirstState] = TmpOperations[i]->EntanglementSpectra[j - TmpOperations[i]->FirstState];
	}
      this->EntanglementSpectrumDimension = TmpOperations[i]->EntanglementSpectrumDimension;
      delete[] TmpOperations[i]->EntanglementSpectra;
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}
  
