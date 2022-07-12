////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of QHE particle hamiltonian precalculation operation        //
//                                                                            //
//                        last modification : 11/03/2003                      //
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
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations

QHEParticlePrecalculationOperation::QHEParticlePrecalculationOperation (AbstractQHEHamiltonian* hamiltonian, bool firstPass)
{
  this->FirstComponent = 0;
  this->Hamiltonian = hamiltonian;
  this->NbrComponent = this->Hamiltonian->GetHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculation;
  this->FirstPass = firstPass;
  this->RequiredMemory = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

QHEParticlePrecalculationOperation::QHEParticlePrecalculationOperation(const QHEParticlePrecalculationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculation;
  this->FirstPass = operation.FirstPass;
  this->RequiredMemory = operation.RequiredMemory;
}
  
// destructor
//

QHEParticlePrecalculationOperation::~QHEParticlePrecalculationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void QHEParticlePrecalculationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* QHEParticlePrecalculationOperation::Clone()
{
  return new QHEParticlePrecalculationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool QHEParticlePrecalculationOperation::RawApplyOperation()
{
  if (this->FirstPass ==  true)
    {
      this->RequiredMemory = this->Hamiltonian->PartialFastMultiplicationMemory(this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->Hamiltonian->PartialEnableFastMultiplication(this->FirstComponent, this->NbrComponent);
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// mpiNodeNbr = provide the additional MPI node ID
// return value = true if no error occurs

bool QHEParticlePrecalculationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr)
{
  long *SegmentIndices=0;
  int TmpNbrThreads = architecture->GetNbrThreads();
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices)==false)
    {
      SegmentIndices = new long[TmpNbrThreads+1];
      int Step = this->NbrComponent / TmpNbrThreads;
      SegmentIndices[0]=this->FirstComponent;
      for (int i=0; i<TmpNbrThreads; ++i)
	SegmentIndices[i]=this->FirstComponent+i*Step;
      SegmentIndices[TmpNbrThreads]=this->FirstComponent+this->NbrComponent;
    }
  QHEParticlePrecalculationOperation** TmpOperations = new QHEParticlePrecalculationOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (QHEParticlePrecalculationOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i+1]-SegmentIndices[i]);
    }
  architecture->SendJobs();
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      if (this->FirstPass ==  true)
	{
	  if (mpiNodeNbr>=0)
	    cout << "node "<<mpiNodeNbr<<" ";
	  cout << "thread "<<i<<" = "<<TmpOperations[i]->RequiredMemory<<endl;
	}
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices)==false)
    delete [] SegmentIndices;
  return true;  
}


