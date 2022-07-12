////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of multiple complex scalar product operation            //
//                                                                            //
//                        last modification : 26/05/2003                      //
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
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/SMPArchitecture.h"

#include <sys/time.h>
#include <stdlib.h>


// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleComplexScalarProductOperation::MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector* rightVectors, 
									     int nbrScalarProduct, Complex* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = rightVectors;
  this->RightVectorsByPointers = 0;
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleComplexScalarProduct;
  this->Strategy = MultipleComplexScalarProductOperation::VectorSubdivision;
  this->ExecutionTime = 0.0;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = array of pointers to the vectors to use for the right hand side of the scalar product
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleComplexScalarProductOperation::MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexVector** rightVectors, 
									     int nbrScalarProduct, Complex* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = rightVectors; 
  this->LeftVector = leftVector;
  this->OperationType = AbstractArchitectureOperation::MultipleComplexScalarProduct;
  this->Strategy = MultipleComplexScalarProductOperation::VectorSubdivision;
  this->ExecutionTime = 0.0;
}

// constructor 
//
// leftVector = pointer to the vector to use for the left hand side of the scalar product
// rightVectors = real matrix where vectors to use for the right hand side of the scalar product are stored
// nbrScalarProduct = number of scalar products that have to be evaluated
// scalarProducts = array where scalar products have to be stored

MultipleComplexScalarProductOperation::MultipleComplexScalarProductOperation(ComplexVector* leftVector, ComplexMatrix& rightVectors, 
									     int nbrScalarProduct, Complex* scalarProducts)
{
  this->FirstComponent = 0;
  this->NbrComponents = leftVector->GetVectorDimension();
  this->NbrScalarProduct = nbrScalarProduct;
  this->ScalarProducts = scalarProducts;
  this->RightVectors = 0; 
  this->RightVectorsByPointers = 0; 
  this->RightVectorMatrix = rightVectors;
  this->LeftVector = leftVector;
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::MultipleComplexScalarProduct;
  this->Strategy = MultipleComplexScalarProductOperation::VectorSubdivision;
}

// copy constructor 
//
// operation = reference on operation to copy

MultipleComplexScalarProductOperation::MultipleComplexScalarProductOperation(const MultipleComplexScalarProductOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponents = operation.NbrComponents;
  this->NbrScalarProduct = operation.NbrScalarProduct;
  this->ScalarProducts = operation.ScalarProducts;
  this->RightVectors = operation.RightVectors;
  this->RightVectorsByPointers = operation.RightVectorsByPointers;
  this->RightVectorMatrix = operation.RightVectorMatrix;
  this->LeftVector = operation.LeftVector;
  this->ExecutionTime = operation.ExecutionTime;
  this->OperationType = AbstractArchitectureOperation::MultipleComplexScalarProduct;
  this->Strategy = operation.Strategy;
}
  
// destructor
//

MultipleComplexScalarProductOperation::~MultipleComplexScalarProductOperation()
{
}
  
// set the array where scalar products have to be stored
//
// scalarProducts = array where scalar products have to be stored

void MultipleComplexScalarProductOperation::SetScalarProducts (Complex* scalarProducts)
{
  this->ScalarProducts = scalarProducts;
}

// set the strategy used to do the scalar products (per vector subdivision or per group)
//
// strategy = flag corresponding to the strategy

void MultipleComplexScalarProductOperation::SetStrategy (int strategy)
{
  this->Strategy = strategy;
}

// set index range of scalar product that have to be calculated
// 
// firstComponent = index of the first component of each partial scalar product for per vector subdivision stategy, 
//                  or index of the first scalar product to evaluate for per group subdivision stategy
// nbrComponent = number of component to take into account for each partial scalar product for per vector subdivision stategy,
//                or number of scalar products that have to be evaluated for per group subdivision stategy

void MultipleComplexScalarProductOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponents = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MultipleComplexScalarProductOperation::Clone()
{
  return new MultipleComplexScalarProductOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool MultipleComplexScalarProductOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);
  if (this->Strategy == MultipleComplexScalarProductOperation::VectorSubdivision)
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
  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MultipleComplexScalarProductOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->LeftVector->GetVectorDimension() / architecture->GetNbrThreads();
  int FirstComponent = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  MultipleComplexScalarProductOperation** TmpOperations = new MultipleComplexScalarProductOperation* [architecture->GetNbrThreads()];
  architecture->SetThreadOperation(this, 0);
  this->SetIndicesRange(FirstComponent, Step);
  FirstComponent += Step;
  for (int i = 1; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (MultipleComplexScalarProductOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(FirstComponent, Step);
      TmpOperations[i]->SetScalarProducts(new Complex [this->GetNbrScalarProduct()]);
      architecture->SetThreadOperation(TmpOperations[i], i);		
      FirstComponent += Step;            
    }
  TmpOperations[ReducedNbrThreads] = (MultipleComplexScalarProductOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(FirstComponent, this->LeftVector->GetVectorDimension() - FirstComponent);  
  TmpOperations[ReducedNbrThreads]->SetScalarProducts(new Complex [this->GetNbrScalarProduct()]);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);		
  architecture->SendJobs();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      sprintf (TmpString, "MultipleComplexScalarProductOperation core operation on SMP id 0 done in %.3f seconds", this->ExecutionTime);
      architecture->AddToLog(TmpString);
      double MinTime = this->ExecutionTime;
      double MaxTime = this->ExecutionTime;
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "MultipleComplexScalarProductOperation core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  if (MinTime > TmpOperations[i]->ExecutionTime)
	    {
	      MinTime = TmpOperations[i]->ExecutionTime;
	    }
	  if (MaxTime < TmpOperations[i]->ExecutionTime)
	    {
	      MaxTime = TmpOperations[i]->ExecutionTime;
	    }
	  architecture->AddToLog(TmpString);
	}
      sprintf (TmpString, "MultipleComplexScalarProductOperation core operation min time=%.3f sec, max time=%.3f sec", MinTime, MaxTime);
      architecture->AddToLog(TmpString);
    }
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
