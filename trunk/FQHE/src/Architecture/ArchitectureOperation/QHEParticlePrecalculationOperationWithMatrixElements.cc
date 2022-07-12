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
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperationWithMatrixElements.h"

#include <sys/time.h>

// debugging switch
// #define DEBUG_QHE_PRECALC_WME

// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations
// tolerance = tolerance for considering matrix elements to be identical

QHEParticlePrecalculationOperationWithMatrixElements::QHEParticlePrecalculationOperationWithMatrixElements (AbstractQHEHamiltonian* hamiltonian, bool firstPass, double tolerance) :
  QHEParticlePrecalculationOperation (hamiltonian, firstPass),
  RealInteractionCoefficients (tolerance),
  ComplexInteractionCoefficients (tolerance)
{
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculationWithMatrixElements;
  this->Tolerance = tolerance;
}

// copy constructor 
//
// operation = reference on operation to copy

QHEParticlePrecalculationOperationWithMatrixElements::QHEParticlePrecalculationOperationWithMatrixElements(const QHEParticlePrecalculationOperationWithMatrixElements& operation):
  QHEParticlePrecalculationOperation(operation),
  RealInteractionCoefficients (operation.Tolerance),
  ComplexInteractionCoefficients (operation.Tolerance)
{
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculationWithMatrixElements;
  this->Tolerance = operation.Tolerance;
}
  
// destructor
//

QHEParticlePrecalculationOperationWithMatrixElements::~QHEParticlePrecalculationOperationWithMatrixElements()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void QHEParticlePrecalculationOperationWithMatrixElements::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  // cout << "QHEParticlePrecalculationOperationWithMatrixElements::SetIndicesRange to "<<firstComponent<<" -> "<< firstComponent + nbrComponent<<endl;
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* QHEParticlePrecalculationOperationWithMatrixElements::Clone()
{
  return new QHEParticlePrecalculationOperationWithMatrixElements (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool QHEParticlePrecalculationOperationWithMatrixElements::RawApplyOperation()
{
  // cout << "RawApplyOperation with "<<  this->FirstComponent << ", "<< this->NbrComponent << endl;

  if (this->FirstPass ==  true)
    {
      this->RequiredMemory = this->Hamiltonian->PartialFastMultiplicationMemory(this->FirstComponent, this->NbrComponent, this->RealInteractionCoefficients, this->ComplexInteractionCoefficients);
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

bool QHEParticlePrecalculationOperationWithMatrixElements::ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr)
{
  // cout << "RawApplyOperation on node "<<mpiNodeNbr<<", firstPass = "<<this->FirstPass<<endl;
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
  QHEParticlePrecalculationOperationWithMatrixElements** TmpOperations = new QHEParticlePrecalculationOperationWithMatrixElements* [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (QHEParticlePrecalculationOperationWithMatrixElements*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i+1]-SegmentIndices[i]);
    }
  architecture->SendJobs();
  
  if (this->FirstPass ==  true) //  merge matrix elements from threads
    {
      this->RealInteractionCoefficients = TmpOperations[0]->RealInteractionCoefficients;
      this->ComplexInteractionCoefficients = TmpOperations[0]->ComplexInteractionCoefficients;
      cout << "Merging matrix elements..."<<endl;
      if (mpiNodeNbr>=0)
	cout << "node "<<mpiNodeNbr<<" ";
      cout << "thread 0 : "<<TmpOperations[0]->ComplexInteractionCoefficients.GetNbrElements()<<" complex, "<<  TmpOperations[0]->RealInteractionCoefficients.GetNbrElements()<<" real"<<endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);
      int StartTimeSecond = TotalStartingTime.tv_sec;      
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  this->RealInteractionCoefficients.MergeArray(TmpOperations[i]->RealInteractionCoefficients);
	  this->ComplexInteractionCoefficients.MergeArray(TmpOperations[i]->ComplexInteractionCoefficients);
	  if (mpiNodeNbr>=0)
	    cout << "node "<<mpiNodeNbr<<" ";
	  cout << "thread "<<i<<" : "<<TmpOperations[i]->ComplexInteractionCoefficients.GetNbrElements()<<" complex, "<<  TmpOperations[i]->RealInteractionCoefficients.GetNbrElements()<<" real"<<endl;
	}
      gettimeofday (&(TotalEndingTime), 0);
      double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);	          
      this->RealInteractionCoefficients.SortEntries();
      this->ComplexInteractionCoefficients.SortEntries();

#ifdef DEBUG_QHE_PRECALC_WME
      // DEBUGGING: testing if all values are present:
      unsigned tmpElementPos;
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
        {
          for (unsigned j=0; j<TmpOperations[i]->RealInteractionCoefficients.GetNbrElements(); ++j)
            if (!this->RealInteractionCoefficients.SearchElement(TmpOperations[i]->RealInteractionCoefficients[j], tmpElementPos))
              cout << "Missing real entry after merging: thread "<<i<<", entry "<<j<<" with value "<< TmpOperations[i]->RealInteractionCoefficients[j] <<endl;
          for (unsigned j=0; j<TmpOperations[i]->ComplexInteractionCoefficients.GetNbrElements(); ++j)
            if (!this->ComplexInteractionCoefficients.SearchElement(TmpOperations[i]->ComplexInteractionCoefficients[j], tmpElementPos))
              cout << "Missing complex entry after merging: thread "<<i<<", entry "<<j<<" with value "<< TmpOperations[i]->RealInteractionCoefficients[j] <<endl;
        }
      this->ComplexInteractionCoefficients.TestAllEntries();
      this->RealInteractionCoefficients.TestAllEntries();
      cout << "Verified all entries"<<endl;
#endif // DEBUG_QHE_PRECALC_WME

      cout << "done merging arrays in "<<Dt<<"s"<<endl;
      if (mpiNodeNbr>=0)
	cout << "node "<<mpiNodeNbr<<" ";
      cout << "merged : "<<this->ComplexInteractionCoefficients.GetNbrElements()<<" complex, "<<  this->RealInteractionCoefficients.GetNbrElements()<<" real"<<endl;
    }
  if (this->FirstPass ==  true)
    {
      cout << "Memory requirements"<<endl;
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  if (mpiNodeNbr>=0)
	    cout << "node "<<mpiNodeNbr<<" ";
	  cout << "thread "<<i<<" = "<<TmpOperations[i]->RequiredMemory<<endl;
	}
    }
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    delete TmpOperations[i];
  delete[] TmpOperations;
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices)==false)
    delete [] SegmentIndices;
  return true;  
}


