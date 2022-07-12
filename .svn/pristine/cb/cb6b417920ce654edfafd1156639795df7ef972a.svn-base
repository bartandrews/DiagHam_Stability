////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of operations that perform the Van der Monde            //
//                      times spinful Slater multiplication                   //
//                                                                            //
//                        last modification : 01/12/2016                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the Hilbert space
// reverseFluxAttachment = use reverse flux attachment
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(BosonOnSphereWithSU2Spin* space, bool reverseFluxAttachment, unsigned long* slaterUp, unsigned long* slaterDown, 
													   double** threeOrbitalOverlaps)
{
  this->Space = (BosonOnSphereWithSU2Spin*) space->Clone();
  this->OutputState = RealVector(this->Space->GetHilbertSpaceDimension(), true);
  this->FirstComponent = 0;
  this->NbrComponent = this->Space->GetNbrParticles();
  this->SlaterUp = new unsigned long[this->Space->NbrBosonsUp];
  for (int i = 0; i < this->Space->NbrBosonsUp; ++i)
    {
      this->SlaterUp[i] = slaterUp[i];
    }
  this->SlaterDown = new unsigned long[this->Space->NbrBosonsDown];
  for (int i = 0; i < this->Space->NbrBosonsDown; ++i)
    {
      this->SlaterDown[i] = slaterDown[i];
    }
  this->Slater2LLUp = 0;
  this->Slater2LLDown = 0;
  this->NbrBosonsLLLUp = this->Space->NbrBosonsUp;
  this->NbrBosonsLLLDown = this->Space->NbrBosonsDown;
  this->ThreeOrbitalOverlaps = new double**[1];
  this->ThreeOrbitalOverlaps[0] = threeOrbitalOverlaps;
  this->ReverseFluxAttachment = reverseFluxAttachment;
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlater;
}

// constructor 
//
// space = pointer to the Hilbert space
// reverseFluxAttachment = use reverse flux attachment
// slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
// slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
// slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
// slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
// nbrBosonsLLLUp = number of spin up bosons in the lowest Landau level
// nbrBosonsLLLDown = number of spin down bosons in the lowest Landau level
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(BosonOnSphereWithSU2Spin* space, bool reverseFluxAttachment, 
													   unsigned long* slaterLLLUp, unsigned long* slater2LLUp, 
													   unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
													   int nbrBosonsLLLUp, int nbrBosonsLLLDown, double*** threeOrbitalOverlaps)
{
  this->Space = (BosonOnSphereWithSU2Spin*) space->Clone();
  this->OutputState = RealVector(this->Space->GetHilbertSpaceDimension(), true);
  this->FirstComponent = 0;
  this->NbrComponent = this->Space->GetNbrParticles();
  this->NbrBosonsLLLUp = nbrBosonsLLLUp;
  this->NbrBosonsLLLDown = nbrBosonsLLLDown;
  this->SlaterUp = new unsigned long[this->NbrBosonsLLLUp];
  for (int i = 0; i < this->NbrBosonsLLLUp; ++i)
    {
      this->SlaterUp[i] = slaterLLLUp[i];
    }
  this->SlaterDown = new unsigned long[this->NbrBosonsLLLDown];
  for (int i = 0; i < this->NbrBosonsLLLDown; ++i)
    {
      this->SlaterDown[i] = slaterLLLDown[i];
    }
  this->Slater2LLUp = new unsigned long[this->Space->NbrBosonsUp - this->NbrBosonsLLLUp];
  for (int i = this->NbrBosonsLLLUp; i < this->Space->NbrBosonsUp; ++i)
    {
      this->Slater2LLUp[i - this->NbrBosonsLLLUp] = slater2LLUp[i - this->NbrBosonsLLLUp];
    }
  this->Slater2LLDown = new unsigned long[this->Space->NbrBosonsDown - this->NbrBosonsLLLDown];
  for (int i = this->NbrBosonsLLLDown; i < this->Space->NbrBosonsDown; ++i)
    {
      this->Slater2LLDown[i - this->NbrBosonsLLLDown] = slater2LLDown[i - this->NbrBosonsLLLDown];
    }
  this->ThreeOrbitalOverlaps = threeOrbitalOverlaps;
  this->ThreeOrbitalOverlaps = threeOrbitalOverlaps;
  this->ReverseFluxAttachment = reverseFluxAttachment;
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlater;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation(const FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Space =  (BosonOnSphereWithSU2Spin*) operation.Space->Clone();
  this->OutputState = operation.OutputState;
  this->ReverseFluxAttachment = operation.ReverseFluxAttachment;
  this->OperationType = AbstractArchitectureOperation::FQHESphereWithSU2SpinVanDerMondeTimesSlater;	
  this->NbrBosonsLLLUp = operation.NbrBosonsLLLUp;
  this->NbrBosonsLLLDown = operation.NbrBosonsLLLDown;
  this->SlaterUp = new unsigned long[this->NbrBosonsLLLUp];
  for (int i = 0; i < this->NbrBosonsLLLUp; ++i)
    {
      this->SlaterUp[i] = operation.SlaterUp[i];
    }
  this->SlaterDown = new unsigned long[this->NbrBosonsLLLDown];
  for (int i = 0; i < this->NbrBosonsLLLDown; ++i)
    {
      this->SlaterDown[i] = operation.SlaterDown[i];
    }
  if (operation.Slater2LLUp != 0)
    {
      this->Slater2LLUp = new unsigned long[this->Space->NbrBosonsUp - this->NbrBosonsLLLUp];
      for (int i = this->NbrBosonsLLLUp; i < this->Space->NbrBosonsUp; ++i)
	{
	  this->Slater2LLUp[i - this->NbrBosonsLLLUp] = operation.Slater2LLUp[i - this->NbrBosonsLLLUp];
	}
      this->Slater2LLDown = new unsigned long[this->Space->NbrBosonsDown - this->NbrBosonsLLLDown];
      for (int i = this->NbrBosonsLLLDown; i < this->Space->NbrBosonsDown; ++i)
	{
	  this->Slater2LLDown[i - this->NbrBosonsLLLDown] = operation.Slater2LLDown[i - this->NbrBosonsLLLDown];
	}
    }
  else
    {
      this->Slater2LLUp = 0;
      this->Slater2LLDown = 0;
    }
  this->ThreeOrbitalOverlaps = operation.ThreeOrbitalOverlaps;
}

// destructor
//

FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::~FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation()
{
  delete this->Space;
  delete[] this->SlaterDown;
  delete[] this->SlaterUp;
  if (this->Slater2LLUp != 0)
    {
      delete[] this->Slater2LLUp;
      delete[] this->Slater2LLDown;
    }
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of components

void FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = (int) firstComponent;
  this->NbrComponent = (int) nbrComponent;
}

// set destination vector 
// 
// destinationVector = new destination vector

void FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::SetDestinationVector (RealVector destinationVector)
{
  this->OutputState = destinationVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::Clone()
{
  return new FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);  
  if (this->Slater2LLUp != 0)
    {
      if (this->ReverseFluxAttachment == true)
	{
	  int LastComponent = this->FirstComponent + this->NbrComponent;
	  for (int j = this->FirstComponent; j < LastComponent; ++j)
	    {
	      this->Space->ReverseVanDerMondeTimesSlater(this->SlaterUp, this->Slater2LLUp, this->SlaterDown, this->Slater2LLDown, this->NbrBosonsLLLUp, this->NbrBosonsLLLDown,
							 this->OutputState, this->ThreeOrbitalOverlaps, j);
	    }
	}
      else
	{
	  cout << "not implemented" << endl;
	}
    }
  else
    {
      if (this->ReverseFluxAttachment == true)
	{
	  int LastComponent = this->FirstComponent + this->NbrComponent;
	  for (int j = this->FirstComponent; j < LastComponent; ++j)
	    {
	      this->Space->ReverseVanDerMondeTimesSlater(this->SlaterUp, this->SlaterDown, this->OutputState, this->ThreeOrbitalOverlaps[0], j);
	    }
	}
      else
	{
	  cout << "not implemented" << endl;
	}
    }
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << "contribution from " << this->FirstComponent << " to " <<  (this->NbrComponent + this->FirstComponent - 1) << " computed in " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation** TmpOperations = new FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation * [architecture->GetNbrThreads()];  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 1; i <  architecture->GetNbrThreads() ; i++)
    {
      RealVector TmpVector (this->OutputState.GetVectorDimension(), true);
      TmpOperations[i]->SetDestinationVector(TmpVector);
    }

  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 

  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      this->OutputState += TmpOperations[i]->OutputState;
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}
