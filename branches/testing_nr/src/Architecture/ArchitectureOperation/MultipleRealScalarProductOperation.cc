////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of multiple real scalar product operation             //
//                                                                            //
//                        last modification : 24/10/2002                      //
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
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/SMPArchitecture.h"


#include <sys/time.h>
#include <stdlib.h>


// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealVector* rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = rightVectors; 
  this->RightVectorsByPointers = 0;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
  this->Strategy = MultipleRealScalarProductOperation::VectorSubdivision;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of pointers to the vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealVector** rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = rightVectors;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
  this->Strategy = MultipleRealScalarProductOperation::VectorSubdivision;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(RealVector* leftVector, RealMatrix& rightVectors, int nbrScalarProduct, double* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = 0;
  this->RightVectorMatrix = rightVectors;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
  this->Strategy = MultipleRealScalarProductOperation::VectorSubdivision;
}

// copy constructor 
//
// operation = reference on operation to copy

MultipleRealScalarProductOperation::MultipleRealScalarProductOperation(const MultipleRealScalarProductOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponents = operation.NbrComponents;
  this->NbrScalarProduct = operation.NbrScalarProduct;
  this->ScalarProducts = operation.ScalarProducts;
  this->RightVectors = operation.RightVectors; 
  this->RightVectorsByPointers = operation.RightVectorsByPointers; 
  this->RightVectorMatrix = operation.RightVectorMatrix;
  this->LeftVector = operation.LeftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleRealScalarProduct;
  this->Strategy = operation.Strategy;
}
  
// destructor
//

MultipleRealScalarProductOperation::~MultipleRealScalarProductOperation()
{
}
  
// set the array where scalar products have to be stored
//
// scalarProducts = array where scalar products have to be stored

void MultipleRealScalarProductOperation::SetScalarProducts (double* scalarProducts)
{
  this->ScalarProducts = scalarProducts;
}

// set the strategy used to do the scalar products (per vector subdivision or per group)
//
// strategy = flag corresponding to the strategy

void MultipleRealScalarProductOperation::SetStrategy (int strategy)
{
  this->Strategy = strategy;
}

// set index range of scalar product that have to be calculated
// 
// firstComponent = index of the first component of each partial scalar product for per vector subdivision stategy, 
//                  or index of the first scalar product to evaluate for per group subdivision stategy
// nbrComponent = number of component to take into account for each partial scalar product for per vector subdivision stategy,
//                or number of scalar products that have to be evaluated for per group subdivision stategy

void MultipleRealScalarProductOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponents = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MultipleRealScalarProductOperation::Clone()
{
  return new MultipleRealScalarProductOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool MultipleRealScalarProductOperation::RawApplyOperation()
{
  if (this->Strategy == MultipleRealScalarProductOperation::VectorSubdivision)
    {
      if (this->RightVectors != 0)
	{
	  for (int i = 0; i < this->NbrScalarProduct; ++i)
	    {
	      this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct(this->RightVectors[i], this->FirstComponent, this->NbrComponents);
	    }
	}
      else
	if (this->RightVectorsByPointers != 0)
	  {
	    for (int i = 0; i < this->NbrScalarProduct; ++i)
	      {
		this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct((*(this->RightVectorsByPointers[i])), this->FirstComponent, this->NbrComponents);
	      }
	  }
	else
	  {
	    for (int i = 0; i < this->NbrScalarProduct; ++i)
	      {
		this->ScalarProducts[i] = this->LeftVector->PartialScalarProduct(this->RightVectorMatrix[i], this->FirstComponent, this->NbrComponents);
	      }
	  }
    }
  else
    {
      int LastScalarProduct = this->FirstComponent + this->NbrComponents;
      if (this->RightVectors != 0)
	{
	  for (int i = this->FirstComponent; i < LastScalarProduct; ++i)
	    {
	      this->ScalarProducts[i] = ((*(this->LeftVector)) * this->RightVectors[i]);
	    }
	}
      else
	if (this->RightVectorsByPointers != 0)
	  {
	    for (int i = this->FirstComponent; i < LastScalarProduct; ++i)
	      {
		this->ScalarProducts[i] = ((*(this->LeftVector)) * (*(this->RightVectorsByPointers[i])));
	      }
	  }
	else
	  {
	    for (int i = this->FirstComponent; i < LastScalarProduct; ++i)
	      {
		this->ScalarProducts[i] = ((*(this->LeftVector)) * this->RightVectorMatrix[i]);
	      }
	  }
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MultipleRealScalarProductOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->LeftVector->GetVectorDimension() / architecture->GetNbrThreads();
  int FirstComponent = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  MultipleRealScalarProductOperation** TmpOperations = new MultipleRealScalarProductOperation* [architecture->GetNbrThreads()];
  architecture->SetThreadOperation(this, 0);
  this->SetIndicesRange(FirstComponent, Step);
  FirstComponent += Step;
  for (int i = 1; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (MultipleRealScalarProductOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(FirstComponent, Step);
      TmpOperations[i]->SetScalarProducts(new double [this->GetNbrScalarProduct()]);
      architecture->SetThreadOperation(TmpOperations[i], i);		
      FirstComponent += Step;            
    }
  TmpOperations[ReducedNbrThreads] = (MultipleRealScalarProductOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(FirstComponent, this->LeftVector->GetVectorDimension() - FirstComponent);  
  TmpOperations[ReducedNbrThreads]->SetScalarProducts(new double [this->GetNbrScalarProduct()]);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);		
  architecture->SendJobs();
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      for (int j = 0; j < this->GetNbrScalarProduct(); ++j)
	this->GetScalarProducts()[j] += TmpOperations[i]->GetScalarProducts()[j];
      delete[] TmpOperations[i]->GetScalarProducts();
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;  
}


// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MultipleRealScalarProductOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
   timeval TotalStartingTime;
   if (architecture->VerboseMode())
     gettimeofday (&TotalStartingTime, 0);
   if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
     this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
   else
     this->RawApplyOperation();
   if (architecture->VerboseMode())
     {
       timeval TotalEndingTime;
       gettimeofday (&TotalEndingTime, 0);
       double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
		     (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
       char TmpString[256];
       sprintf (TmpString, "MultipleRealScalarProductOperation core operation done in %.3f seconds", Dt);
       architecture->AddToLog(TmpString, true);
     }
   return true;
#else
   return false;
#endif
}
