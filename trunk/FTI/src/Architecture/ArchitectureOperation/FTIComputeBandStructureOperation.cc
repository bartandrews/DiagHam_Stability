////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FTI band structure calculation operation             //
//                                                                            //
//                        last modification : 26/09/2012                      //
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
#include "Architecture/ArchitectureOperation/FTIComputeBandStructureOperation.h"

#include <sys/time.h>


// constructor 
//
// tightBindingModel = pointer to the tight binding model

FTIComputeBandStructureOperation::FTIComputeBandStructureOperation (AbstractTightBindingModel* tightBindingModel)
{
  this->TightBindingModel = tightBindingModel;
  this->OperationType = AbstractArchitectureOperation::FTIComputeBandStructureOperation;
  this->FirstComponent = 0;
  this->NbrComponent = this->TightBindingModel->GetNbrStatePerBand();
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = this->TightBindingModel->GetNbrStatePerBand();
}

// copy constructor 
//
// operation = reference on operation to copy

FTIComputeBandStructureOperation::FTIComputeBandStructureOperation(const FTIComputeBandStructureOperation& operation)
{
  this->OperationType = AbstractArchitectureOperation::FTIComputeBandStructureOperation;
  this->TightBindingModel= operation.TightBindingModel;
}
  
// destructor
//

FTIComputeBandStructureOperation::~FTIComputeBandStructureOperation()
{
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FTIComputeBandStructureOperation::Clone()
{
  return new FTIComputeBandStructureOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FTIComputeBandStructureOperation::RawApplyOperation()
{
  if (this->LargeNbrComponent > 0l)
    {
      this->TightBindingModel->CoreComputeBandStructure(this->LargeFirstComponent, this->LargeNbrComponent);
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FTIComputeBandStructureOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long TmpNbrThreads = ((long) architecture->GetNbrThreads());
  if(this->LargeNbrComponent < TmpNbrThreads)
    TmpNbrThreads = this->LargeNbrComponent;
  long Step = this->LargeNbrComponent / ((long) TmpNbrThreads);
  FTIComputeBandStructureOperation** TmpOperations = new FTIComputeBandStructureOperation* [TmpNbrThreads];
  long DecTmpNbrThreads = TmpNbrThreads - 1;
  for (int i = 0; i < DecTmpNbrThreads ; i++)
    {
      TmpOperations[i] = (FTIComputeBandStructureOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(i * Step, Step);
    }
  TmpOperations[DecTmpNbrThreads] = (FTIComputeBandStructureOperation *) this->Clone();
  architecture->SetThreadOperation(TmpOperations[DecTmpNbrThreads], DecTmpNbrThreads);
  TmpOperations[DecTmpNbrThreads]->SetIndicesRange(DecTmpNbrThreads * Step, 
						   this->LargeNbrComponent - DecTmpNbrThreads * Step);
  
  architecture->SendJobs(TmpNbrThreads);
  
  for (int i = 0; i <  TmpNbrThreads; ++i)
    delete TmpOperations[i];
  delete[] TmpOperations;
  
  return true;
}
  
bool FTIComputeBandStructureOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
  return this->RawApplyOperation();
}
