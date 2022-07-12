////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of abstract hamiltonian precalculation operation        //
//                                                                            //
//                        last modification : 13/11/2003                      //
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
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Architecture/SMPArchitecture.h"


// destructor
//

AbstractPrecalculationOperation::~AbstractPrecalculationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AbstractPrecalculationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
  this->LargeFirstComponent = (long) firstComponent;
  this->LargeNbrComponent = (long) nbrComponent;
}

// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AbstractPrecalculationOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->LargeFirstComponent = firstComponent;
  this->LargeNbrComponent = nbrComponent;
  if (this->LargeFirstComponent < (1l << 30))
    this->FirstComponent = (int) this->LargeFirstComponent;    
  else
    this->FirstComponent = 0;
  if (this->LargeNbrComponent < (1l << 30))
    this->NbrComponent = (int) this->LargeNbrComponent;    
  else
    this->NbrComponent = 0;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// mpiNodeNbr = provide the additional MPI node ID
// return value = true if no error occurs

bool AbstractPrecalculationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  AbstractPrecalculationOperation** TmpOperations = new AbstractPrecalculationOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (AbstractPrecalculationOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (AbstractPrecalculationOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AbstractPrecalculationOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  this->SetIndicesRange(TmpMinimumIndex, (TmpMaximumIndex - TmpMinimumIndex + 1));
  if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
    this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture(), architecture->GetNodeNbr());
  else
    this->RawApplyOperation();
  return true;
}
  
// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AbstractPrecalculationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture) 
{
  return this->ArchitectureDependentApplyOperation(architecture, -1);
}

