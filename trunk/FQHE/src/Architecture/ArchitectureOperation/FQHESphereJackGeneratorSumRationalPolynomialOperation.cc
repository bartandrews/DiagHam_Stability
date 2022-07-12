////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FQHE on disk quasihole propagator operation          //
//                                                                            //
//                        last modification : 08/03/2009                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorSumRationalPolynomialOperation.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"


// constructor 
//
// space = pointer to the Hilbert space to use
// invAlpha = inverse of the Jack polynomial alpha coefficient
// rootPartition = root partition (in fermionic binary representation)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// rhoArray = rho factor associated to each state
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// fermionicFlag = true if we are dealing with fermions
// symmetricFlag = true if the state is Lz<->-Lz symmetric

FQHESphereJackGeneratorSumRationalPolynomialOperation::FQHESphereJackGeneratorSumRationalPolynomialOperation (long index, LongRationalPolynomial* numerators, LongRationalPolynomial* denominators, long* connectedIndices, long* connectedCoefficients, int nbrConnected)
{
  this->Numerators = numerators;
  this->Denominators = denominators;
  this->ConnectedIndices = connectedIndices;
  this->ConnectedCoefficients = connectedCoefficients;
  this->LocalNbrRationalPolynomials = nbrConnected; 
  this->LocalNumerators = new LongRationalPolynomial*[(this->LocalNbrRationalPolynomials >> 1) + 1];
  this->LocalDenominators = new LongRationalPolynomial*[(this->LocalNbrRationalPolynomials >> 1) + 1];
  this->LocalNumerators[0] = &(this->Numerators[index]);
  this->LocalDenominators[0] = &(this->Denominators[index]);
  this->LocalShift = 0;
  this->LocalIndex = 0;
  this->LocalStep = 1;
  this->LocalNbrJobs = 1;
  this->FirstPassFlag = true;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJackGenerator;
  this->Flag.Initialize();
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereJackGeneratorSumRationalPolynomialOperation::FQHESphereJackGeneratorSumRationalPolynomialOperation(const FQHESphereJackGeneratorSumRationalPolynomialOperation& operation)
{
  this->Numerators = operation.Numerators;
  this->Denominators = operation.Denominators;
  this->ConnectedIndices = operation.ConnectedIndices;
  this->ConnectedCoefficients = operation.ConnectedCoefficients;
  this->LocalNbrRationalPolynomials = operation.LocalNbrRationalPolynomials; 
  this->LocalNumerators = operation.LocalNumerators;
  this->LocalDenominators = operation.LocalDenominators;
  this->LocalShift = operation.LocalShift;
  this->LocalIndex = operation.LocalIndex;
  this->LocalStep = operation.LocalStep;
  this->LocalNbrJobs = operation.LocalNbrJobs;
  this->FirstPassFlag = operation.FirstPassFlag;
  this->OperationType = AbstractArchitectureOperation::FQHESphereJackGenerator;
  this->Flag = operation.Flag;
}
  
// destructor
//

FQHESphereJackGeneratorSumRationalPolynomialOperation::~FQHESphereJackGeneratorSumRationalPolynomialOperation()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 1; i < this->LocalNbrJobs; ++i)
	{
	  delete this->LocalNumerators[i];
	  delete this->LocalDenominators[i];
	}
      delete[] this->LocalNumerators;
      delete[] this->LocalDenominators;
    }
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereJackGeneratorSumRationalPolynomialOperation::Clone()
{
  return new FQHESphereJackGeneratorSumRationalPolynomialOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereJackGeneratorSumRationalPolynomialOperation::RawApplyOperation()
{
  if (this->FirstPassFlag == true)
    {
      LongRationalPolynomial* TmpLocalNumerator = this->LocalNumerators[this->LocalIndex];
      LongRationalPolynomial* TmpLocalDenominator = this->LocalDenominators[this->LocalIndex];
      (*TmpLocalNumerator) = this->Numerators[this->ConnectedIndices[this->LocalShift]];
      (*TmpLocalDenominator) = this->Denominators[this->ConnectedIndices[this->LocalShift]];
      (*TmpLocalNumerator) *= this->ConnectedCoefficients[this->LocalShift];
      if (this->LocalNbrRationalPolynomials > 1)
	{
	  for (int i = 1; i < this->LocalNbrRationalPolynomials; ++i)
	    {
	      LongRationalPolynomial TmpNumerator (this->Numerators[this->ConnectedIndices[i + this->LocalShift]], (*TmpLocalDenominator));
	      TmpNumerator *= this->ConnectedCoefficients[i + this->LocalShift];
	      LongRationalPolynomial& TmpPolynomial = this->Denominators[this->ConnectedIndices[i + this->LocalShift]];
	      (*TmpLocalNumerator) *= TmpPolynomial;
	      (*TmpLocalNumerator) += TmpNumerator;
	      LongRationalPolynomial TmpDenominator (0);
	      TmpDenominator[0] = 1l;
	      for (int j = 0; j < TmpLocalDenominator->GetPolynomialDegree(); ++j)
		{
		  TmpLocalNumerator->PolynomialEvaluate(TmpLocalDenominator->PolynomialRationalRoot(j), this->TmpRational);	
		  if (this->TmpRational.IsZero())
		    {
		      TmpLocalNumerator->LocalMonomialDivision(TmpLocalDenominator->PolynomialRationalRoot(j));
		    }
		  else
		    {
		      TmpDenominator.LocalMonomialMultiplication(TmpLocalDenominator->PolynomialRationalRoot(j));
		    }
		}
	      (*TmpLocalDenominator) = TmpDenominator;
	      for (int j = 0; j < TmpPolynomial.GetPolynomialDegree(); ++j)
		{
		  TmpLocalNumerator->PolynomialEvaluate(TmpPolynomial.PolynomialRationalRoot(j), this->TmpRational);
		  if (this->TmpRational.IsZero())
		    {
		      TmpLocalNumerator->LocalMonomialDivision(TmpPolynomial.PolynomialRationalRoot(j));
		    }
		  else
		    {
		      TmpLocalDenominator->LocalMonomialMultiplication(TmpPolynomial.PolynomialRationalRoot(j));
		    }
		}
	    }
	}
    }
  else
    {
      if (this->LocalNbrRationalPolynomials > 1)
	{
	  LongRationalPolynomial& TmpLocalNumerator = *(this->LocalNumerators[this->LocalIndex]);
	  LongRationalPolynomial& TmpLocalDenominator = *(this->LocalDenominators[this->LocalIndex]);
	  for (int i = 1; i < this->LocalNbrRationalPolynomials; ++i)
	    {
	      LongRationalPolynomial TmpNumerator (*(this->LocalNumerators[(i * this->LocalStep) + this->LocalIndex]), TmpLocalDenominator);
	      LongRationalPolynomial& TmpPolynomial = *(this->LocalDenominators[(i * this->LocalStep) + this->LocalIndex]);
	      TmpLocalNumerator *= TmpPolynomial;
	      TmpLocalNumerator += TmpNumerator;
	      LongRationalPolynomial TmpDenominator (0);
	      TmpDenominator[0] = 1l;
	      for (int j = 0; j < TmpLocalDenominator.GetPolynomialDegree(); ++j)
		{
		  TmpLocalNumerator.PolynomialEvaluate(TmpLocalDenominator.PolynomialRationalRoot(j), this->TmpRational);	
		  if (this->TmpRational.IsZero())
		    {
		      TmpLocalNumerator.LocalMonomialDivision(TmpLocalDenominator.PolynomialRationalRoot(j));
		    }
		  else
		    {
		      TmpDenominator.LocalMonomialMultiplication(TmpLocalDenominator.PolynomialRationalRoot(j));
		    }
		}
	      TmpLocalDenominator = TmpDenominator;
	      for (int j = 0; j < TmpPolynomial.GetPolynomialDegree(); ++j)
		{
		  TmpLocalNumerator.PolynomialEvaluate(TmpPolynomial.PolynomialRationalRoot(j), this->TmpRational);
		  if (this->TmpRational.IsZero())
		    {
		      TmpLocalNumerator.LocalMonomialDivision(TmpPolynomial.PolynomialRationalRoot(j));
		    }
		  else
		    {
		      TmpLocalDenominator.LocalMonomialMultiplication(TmpPolynomial.PolynomialRationalRoot(j));
		    }
		}
	    }  
	}
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereJackGeneratorSumRationalPolynomialOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  if (this->LocalNbrRationalPolynomials <= 3)
    {
      this->RawApplyOperation();
      return true;
    }
  this->LocalNbrJobs = architecture->GetNbrThreads();
  while ((this->LocalNbrRationalPolynomials / this->LocalNbrJobs) <= 1)
    --this->LocalNbrJobs;
  for (int i = 1; i < this->LocalNbrJobs; ++i)
    {
      this->LocalNumerators[i] = new LongRationalPolynomial;
      this->LocalDenominators[i] = new LongRationalPolynomial;
    }
  int Step = this->LocalNbrRationalPolynomials / this->LocalNbrJobs;
  
  FQHESphereJackGeneratorSumRationalPolynomialOperation** LocalOperations = new FQHESphereJackGeneratorSumRationalPolynomialOperation* [this->LocalNbrJobs];
  for (int i = 0; i < this->LocalNbrJobs; ++i)
    LocalOperations[i] = (FQHESphereJackGeneratorSumRationalPolynomialOperation*) this->Clone();
  
  int ReducedNbrThreads = this->LocalNbrJobs - 1;
  int TmpFirstComponent = 0;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      LocalOperations[i]->LocalShift = TmpFirstComponent;
      LocalOperations[i]->LocalNbrRationalPolynomials = Step;
      LocalOperations[i]->LocalIndex = i;
      architecture->SetThreadOperation(LocalOperations[i], i);
      TmpFirstComponent += Step;
    }
  LocalOperations[ReducedNbrThreads]->LocalShift = TmpFirstComponent;
  LocalOperations[ReducedNbrThreads]->LocalNbrRationalPolynomials = this->LocalNbrRationalPolynomials - TmpFirstComponent;
  LocalOperations[ReducedNbrThreads]->LocalIndex = ReducedNbrThreads;
  architecture->SetThreadOperation(LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs(this->LocalNbrJobs);
  
  for (int i = 0; i < this->LocalNbrJobs; ++i)
    LocalOperations[i]->FirstPassFlag = false;
    
  int TmpNbrJobs = this->LocalNbrJobs;
  Step = 1;
  while ((TmpNbrJobs >> 1) > 2)
    {
      ReducedNbrThreads = (TmpNbrJobs >> 1) - 1;
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  LocalOperations[i]->LocalNbrRationalPolynomials = 2;
	  LocalOperations[i]->LocalIndex = i * (Step << 1);
	  LocalOperations[i]->LocalStep = Step;
	}
      LocalOperations[ReducedNbrThreads]->LocalNbrRationalPolynomials = 2 + (TmpNbrJobs & 1);
      LocalOperations[ReducedNbrThreads]->LocalIndex = ReducedNbrThreads * (Step << 1);
      LocalOperations[ReducedNbrThreads]->LocalStep = Step;
      TmpNbrJobs >>= 1;
      architecture->SendJobs(TmpNbrJobs);
      Step <<= 1;
    }
  LocalOperations[0]->LocalNbrRationalPolynomials = TmpNbrJobs;
  LocalOperations[0]->LocalIndex = 0;
  LocalOperations[0]->LocalStep = Step;
  LocalOperations[0]->RawApplyOperation();

  for (int i = 0; i < this->LocalNbrJobs; ++i)
    delete LocalOperations[i];
  delete[] LocalOperations;
  return true;
}
  
