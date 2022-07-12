////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract scalar sum  operation                  //
//                                                                            //
//                        last modification : 29/07/2003                      //
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
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "Architecture/SMPArchitecture.h"


// default constructor
//

AbstractScalarSumOperation::AbstractScalarSumOperation()
{
  this->LongRationalScalars = 0;
}
  
// destructor
//

AbstractScalarSumOperation::~AbstractScalarSumOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AbstractScalarSumOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
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

void AbstractScalarSumOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
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
// return value = true if no error occurs

bool AbstractScalarSumOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->GetLargeDimension() / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = 0l;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  AbstractScalarSumOperation** TmpOperations = new AbstractScalarSumOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (AbstractScalarSumOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (AbstractScalarSumOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->GetLargeDimension() - TmpFirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  if (this->NbrScalars > 1)
    {
      if (this->Scalars != 0)
	{
	  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	    {
	      for (int j = 0; j < this->NbrScalars; ++j)
		this->GetScalar(j) += TmpOperations[i]->GetScalar(j);
	      delete TmpOperations[i];
	    }
	}
      else
	{
	  if (this->LongRationalScalars != 0)
	    {
	      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
		{
		  for (int j = 0; j < this->NbrScalars; ++j)
		    this->GetLongRationalScalar(j) += TmpOperations[i]->GetLongRationalScalar(j);
		  delete TmpOperations[i];
		}
	    }
	}
    }
  else
    {
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  this->GetScalar() += TmpOperations[i]->GetScalar();
	  this->GetLongRationalScalar() += TmpOperations[i]->GetLongRationalScalar();
	  delete TmpOperations[i];
	}
    }
  delete[] TmpOperations;
  return true;
}
