////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of Abstract spin chain with translations              //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"

#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// destructor
//

AbstractDoubledSpinChainWithTranslations::~AbstractDoubledSpinChainWithTranslations () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if(this->RescalingFactors != 0 )
	{ 
	  for (int i = 1; i <= this->MaxXMomentum; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	  delete[] this->CompatibilityWithMomentum; 
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

AbstractDoubledSpinChainWithTranslations & AbstractDoubledSpinChainWithTranslations::operator = (const AbstractDoubledSpinChainWithTranslations & chain)
{
//  AbstractDoubledSpinChain::operator = (chain);

  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->CompatibilityWithMomentum;
	  for (int i = 1; i <= this->ChainLength; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	}
    }  
  if (chain.ChainLength != 0)
    {
  
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
   }
  else
    {
      this->Momentum = 0;
      this->MaxXMomentum = 0;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// create precalculation tables
//

void AbstractDoubledSpinChainWithTranslations::CreatePrecalculationTable()
{
  int TmpPeriodicity = this->MaxXMomentum ;
  this->CompatibilityWithMomentum = new bool [TmpPeriodicity + 1];
  for (int i = 0; i <= TmpPeriodicity; ++i)
    if (((i * this->Momentum) % TmpPeriodicity) == 0)
      this->CompatibilityWithMomentum[i] = true;
    else
      this->CompatibilityWithMomentum[i] = false;

  this->RescalingFactors = new double* [TmpPeriodicity + 1];
  for (int i = 1; i <= TmpPeriodicity; ++i)
    {
      this->RescalingFactors[i] = new double [TmpPeriodicity + 1];
      for (int j = 1; j <= TmpPeriodicity; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}
