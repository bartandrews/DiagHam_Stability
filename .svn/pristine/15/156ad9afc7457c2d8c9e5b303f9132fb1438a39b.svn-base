////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                    for torus with magnetic translations                    //
//                                                                            //
//                        last modification : 13/06/2011                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"

#include <sys/time.h>


// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithMagneticTranslations* fullSpace, ParticleOnTorusWithMagneticTranslations* destinationSpace, ParticleOnTorusWithMagneticTranslations* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnTorusWithMagneticTranslations*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) complementarySpace->Clone();
  this->SpinfulFullSpace = 0;
  this->SpinfulDestinationHilbertSpace = 0; 
  this->SpinfulComplementaryHilbertSpace = 0;
  this->ComplexGroundState = groundState;
  this->ComplexDensityMatrix = densityMatrix;
  this->ComplexGroundStates = 0;
  this->GroundStateWeights = 0;
  this->NbrGroundStates = 0;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
//   cout << "subsystem dimension = " << this->DestinationHilbertSpace->GetLargeHilbertSpaceDimension() << " "  << this->DestinationHilbertSpace->GetHilbertSpaceDimension() << endl; 
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// constructor for the spinful case
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithSpinAndMagneticTranslations* fullSpace, ParticleOnTorusWithSpinAndMagneticTranslations* destinationSpace, ParticleOnTorusWithSpinAndMagneticTranslations* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix)
{
  this->SpinfulFullSpace  = (ParticleOnTorusWithSpinAndMagneticTranslations*) fullSpace->Clone();
  this->SpinfulDestinationHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) destinationSpace->Clone();
  this->SpinfulComplementaryHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) complementarySpace->Clone();
  this->FullSpace = 0;
  this->DestinationHilbertSpace = 0; 
  this->ComplementaryHilbertSpace = 0;
  this->ComplexGroundState = groundState;
  this->ComplexGroundStates = 0;
  this->GroundStateWeights = 0;
  this->NbrGroundStates = 0;
  this->ComplexDensityMatrix = densityMatrix;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->SpinfulComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
//   cout << "subsystem dimension = " << this->DestinationHilbertSpace->GetLargeHilbertSpaceDimension() << " "  << this->DestinationHilbertSpace->GetHilbertSpaceDimension() << endl; 
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// constructor when using a sum of projectors
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithMagneticTranslations* fullSpace, 
											       ParticleOnTorusWithMagneticTranslations* destinationSpace, 
											       ParticleOnTorusWithMagneticTranslations* complementarySpace, 
											       int nbrGroundStates, ComplexVector* groundStates, double* weights, 
											       HermitianMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnTorusWithMagneticTranslations*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) complementarySpace->Clone();
  this->SpinfulFullSpace = 0;
  this->SpinfulDestinationHilbertSpace = 0; 
  this->SpinfulComplementaryHilbertSpace = 0;
  this->ComplexDensityMatrix = densityMatrix;
  this->NbrGroundStates = nbrGroundStates;
  this->ComplexGroundStates = new ComplexVector[this->NbrGroundStates];
  this->GroundStateWeights = new double[this->NbrGroundStates];
  for (int i = 0; i < this->NbrGroundStates; ++i)
    {
      this->ComplexGroundStates[i] = groundStates[i];
      this->GroundStateWeights[i] = weights[i];
    }
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
//   cout << "subsystem dimension = " << this->DestinationHilbertSpace->GetLargeHilbertSpaceDimension() << " "  << this->DestinationHilbertSpace->GetHilbertSpaceDimension() << endl; 
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// constructor for the spinful case when using a sum of projectors
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithSpinAndMagneticTranslations* fullSpace, 
											       ParticleOnTorusWithSpinAndMagneticTranslations* destinationSpace, 
											       ParticleOnTorusWithSpinAndMagneticTranslations* complementarySpace, 
											       int nbrGroundStates, ComplexVector* groundStates, double* weights, 
											       HermitianMatrix& densityMatrix)
{
  this->SpinfulFullSpace  = (ParticleOnTorusWithSpinAndMagneticTranslations*) fullSpace->Clone();
  this->SpinfulDestinationHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) destinationSpace->Clone();
  this->SpinfulComplementaryHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) complementarySpace->Clone();
  this->FullSpace = 0;
  this->DestinationHilbertSpace = 0; 
  this->ComplementaryHilbertSpace = 0;
  this->NbrGroundStates = nbrGroundStates;
  this->ComplexGroundStates = new ComplexVector[this->NbrGroundStates];
  this->GroundStateWeights = new double[this->NbrGroundStates];
  for (int i = 0; i < this->NbrGroundStates; ++i)
    {
      this->ComplexGroundStates[i] = groundStates[i];
      this->GroundStateWeights[i] = weights[i];
    }
  this->ComplexDensityMatrix = densityMatrix;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->SpinfulComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
//   cout << "subsystem dimension = " << this->DestinationHilbertSpace->GetLargeHilbertSpaceDimension() << " "  << this->DestinationHilbertSpace->GetHilbertSpaceDimension() << endl; 
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(const FQHETorusParticleEntanglementSpectrumOperation& operation)
{
  this->NbrGroundStates = operation.NbrGroundStates;
  if (this->NbrGroundStates > 0)
    {
      this->ComplexGroundStates = new ComplexVector[this->NbrGroundStates];
      this->GroundStateWeights = new double[this->NbrGroundStates];
      for (int i = 0; i < this->NbrGroundStates; ++i)
	{
	  this->ComplexGroundStates[i] = operation.ComplexGroundStates[i];
	  this->GroundStateWeights[i] = operation.GroundStateWeights[i];
	}
    }
  if (operation.SpinfulFullSpace == 0)
    {
      this->FullSpace  = (ParticleOnTorusWithMagneticTranslations*) operation.FullSpace->Clone();
      this->DestinationHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) operation.DestinationHilbertSpace->Clone();
      this->ComplementaryHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) operation.ComplementaryHilbertSpace->Clone();
      this->SpinfulFullSpace = 0;
      this->SpinfulDestinationHilbertSpace = 0; 
      this->SpinfulComplementaryHilbertSpace = 0;
    }
  else
    {
      this->FullSpace = 0;
      this->DestinationHilbertSpace = 0; 
      this->ComplementaryHilbertSpace = 0;
      this->SpinfulFullSpace  = (ParticleOnTorusWithSpinAndMagneticTranslations*) operation.SpinfulFullSpace->Clone();
      this->SpinfulDestinationHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) operation.SpinfulDestinationHilbertSpace->Clone();
      this->SpinfulComplementaryHilbertSpace = (ParticleOnTorusWithSpinAndMagneticTranslations*) operation.SpinfulComplementaryHilbertSpace->Clone();
    }
  this->ComplexGroundState = operation.ComplexGroundState;
  this->ComplexDensityMatrix = operation.ComplexDensityMatrix;
  this->NbrNonZeroElements = operation.NbrNonZeroElements;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}
  
// destructor
//

FQHETorusParticleEntanglementSpectrumOperation::~FQHETorusParticleEntanglementSpectrumOperation()
{
  if (this->LocalOperations != 0)
    {
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	{
	  delete this->LocalOperations[i];
	}
      delete[] this->LocalOperations;
    }
  if (this->NbrGroundStates > 0)
    {
      delete[] this->ComplexGroundStates;
      delete[] this->GroundStateWeights;
    }
  if (this->SpinfulFullSpace == 0)
    {
      delete this->FullSpace;
      delete this->DestinationHilbertSpace;
      delete this->ComplementaryHilbertSpace;
    }
  else
    {
      delete this->SpinfulFullSpace;
      delete this->SpinfulDestinationHilbertSpace;
      delete this->SpinfulComplementaryHilbertSpace;
    }
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHETorusParticleEntanglementSpectrumOperation::Clone()
{
  return new FQHETorusParticleEntanglementSpectrumOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusParticleEntanglementSpectrumOperation::RawApplyOperation()
{
  if (this->NbrGroundStates == 0)
    {
      if (this->SpinfulFullSpace == 0)
	{
	  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->ComplexGroundState, &this->ComplexDensityMatrix);
	}
      else
	{
	  this->NbrNonZeroElements = this->SpinfulFullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->SpinfulComplementaryHilbertSpace, this->SpinfulDestinationHilbertSpace, this->ComplexGroundState, &this->ComplexDensityMatrix);
	}
    }
  else
    {
      if (this->SpinfulFullSpace == 0)
	{
	  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->NbrGroundStates, this->ComplexGroundStates, this->GroundStateWeights, &this->ComplexDensityMatrix);
	}
      else
	{
	  this->NbrNonZeroElements = this->SpinfulFullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->SpinfulComplementaryHilbertSpace, this->SpinfulDestinationHilbertSpace, this->NbrGroundStates, this->ComplexGroundStates, this->GroundStateWeights, &this->ComplexDensityMatrix);
	}
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusParticleEntanglementSpectrumOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  if (this->LocalOperations == 0)
    {
      this->NbrLocalOperations = architecture->GetNbrThreads();
      this->LocalOperations = new FQHETorusParticleEntanglementSpectrumOperation* [this->NbrLocalOperations];
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	this->LocalOperations[i] = (FQHETorusParticleEntanglementSpectrumOperation*) this->Clone();
    }
  int ReducedNbrThreads = this->NbrLocalOperations - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      if (this->SpinfulComplementaryHilbertSpace == 0)
	this->LocalOperations[i]->ComplexDensityMatrix = HermitianMatrix(this->DestinationHilbertSpace->GetHilbertSpaceDimension(), true);
      else
	this->LocalOperations[i]->ComplexDensityMatrix = HermitianMatrix(this->SpinfulDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
      architecture->SetThreadOperation(this->LocalOperations[i], i);
      TmpFirstComponent += Step;
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[ReducedNbrThreads]->ComplexDensityMatrix += this->LocalOperations[i]->ComplexDensityMatrix;
      this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements += this->LocalOperations[i]->NbrNonZeroElements;
    }
  this->NbrNonZeroElements = this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements;
  return true;
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusParticleEntanglementSpectrumOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->LargeFirstComponent = TmpMinimumIndex;  
  this->LargeNbrComponent = (TmpMaximumIndex - TmpMinimumIndex + 1l);  
  this->FirstComponent = (int) TmpMinimumIndex;  
  this->NbrComponent = (int) (TmpMaximumIndex - TmpMinimumIndex + 1l);
  timeval TotalStartingTime;
  if (architecture->VerboseMode())
    gettimeofday (&TotalStartingTime, 0);
  if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
    {
      this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
    }
  else
    this->RawApplyOperation();
  if (architecture->VerboseMode())
    {
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
      char TmpString[256];
      sprintf (TmpString, "FQHETorusParticleEntanglementSpectrumOperation core operation done in %.3f seconds", Dt);
      architecture->AddToLog(TmpString);
    }
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    gettimeofday (&TotalStartingTime, 0);
  architecture->SumMatrix(this->ComplexDensityMatrix);
//  architecture->BroadcastMatrix(this->ComplexDensityMatrix);
  if (architecture->IsMasterNode() == false)
    {
      this->ComplexDensityMatrix = HermitianMatrix(2, true);
      this->ComplexDensityMatrix.SetMatrixElement(0, 0, 0.5);
      this->ComplexDensityMatrix.SetMatrixElement(1, 1, 0.5);
    }
  if (architecture->IsMasterNode())
    {
      long TmpValues;;
      int TmpNbrValues = 1 ;
      for (int i = 0; i < architecture->GetNbrSlaveNodes(); ++i)
	{
	  architecture->ReceiveFromSlave(i, &TmpValues, TmpNbrValues);
	  this->NbrNonZeroElements += TmpValues;
	}
      TmpValues = this->NbrNonZeroElements;      
      architecture->BroadcastToSlaves(&TmpValues, TmpNbrValues);
    }
  else
    {
      long TmpValues = this->NbrNonZeroElements;
      int TmpNbrValues = 1 ;     
      architecture->SendToMaster(&TmpValues, TmpNbrValues);
      architecture->BroadcastToSlaves(&TmpValues, TmpNbrValues);
      this->NbrNonZeroElements = TmpValues;
   }
  if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
    {
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
      char TmpString[256];
      sprintf (TmpString, "FQHETorusParticleEntanglementSpectrumOperation sum operation done in %.3f seconds", Dt);
      architecture->AddToLog(TmpString, true);
    }
  return true;
#else
  return this->RawApplyOperation();
#endif
}
