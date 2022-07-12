////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                     class author: Cecile Repellin                          //
//                                                                            //
//                                                                            //
//              class of fermions on a honeycomb lattice in real space        //
//                        with plaquette cluster exclusion                    //
//                                                                            //
//                        last modification : 01/04/2018                      //
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
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion ()
{
  this->NbrSitesX = 0;
  this->NbrSitesY = 0;
  this->ListIndicesPerPlaquette = 0;
  this->LargerIndicesInPlaquette = 0;
  this->NbrLargerIndicesInPlaquette = 0;
  this->SzFlag = false;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// memory = amount of memory granted for precalculations

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion (int nbrFermions, int nbrSitesX, int nbrSitesY, unsigned long memory)
{  
  this->SzFlag = false;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSitesX = nbrSitesX;
  this->NbrSitesY = nbrSitesY;  
  this->NbrSite = 2 * this->NbrSitesX * this->NbrSitesY;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  
  
  this->InitializeHexagonArrays();
   
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSite - 1, 0x0ul);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  unsigned long TmpUnfilteredDimension = (FermionOnLatticeWithSpinRealSpace::EvaluateHilbertSpaceDimension(this->NbrFermions));
  cout << "Hilbert space dimension unfiltered = " << TmpUnfilteredDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, 0x0ul, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
         
      
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  if (this->FindStateIndex(this->StateDescription[i],  this->StateHighestBit[i]) != i)
	    {
	      cout << i << " " << this->FindStateIndex(this->StateDescription[i],  this->StateHighestBit[i]) << endl;
	    }
	}


#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}


// constructor for preserved Sz
// 
// nbrFermions = number of fermions
// totalSpin = twice the total spin value
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// memory = amount of memory granted for precalculations

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion (int nbrFermions, int totalSpin, int nbrSitesX, int nbrSitesY, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSitesX = nbrSitesX;
  this->NbrSitesY = nbrSitesY;  
  this->NbrSite = 2 * this->NbrSitesX * this->NbrSitesY;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  
  
  this->InitializeHexagonArrays();
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp, this->NbrSite - 1, 0x0ul);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  int TmpUnfilteredDimension = (FermionOnLatticeWithSpinRealSpace::EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp));
  cout << "Hilbert space dimension unfiltered = " << TmpUnfilteredDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrFermionsUp, this->NbrSite - 1, 0x0ul, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
      
//       for (int i = 0; i < this->NbrSite; ++i)
//       {
// 	for (int j = 0; j  <this->NbrSite; ++j)
// 	{
// 	  double Tmpcoef = this->AuAd(453, i, j);
// 	  double Tmpcoef1;
// 	  int Index = this->AduAdd(i, j, Tmpcoef1);
// 	  if (Index == 453)
// 	    cout << i << " " << j << " : " << Tmpcoef << " " << Tmpcoef1 << " " << Index << endl;
// 	}
//       }
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//       {
// 	  cout << i << " : " ;
// 	  this->PrintState(cout, i);
// 	  cout << endl;
//       }
//       
//       this->StateDescription = new unsigned long [TmpUnfilteredDimension];
//       FermionOnLatticeWithSpinRealSpace::GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0l); 
// //       test the consistency of the constraint
//       int TmpHilbertSpaceDimension = 0;
//       for (int i = 0; i < TmpUnfilteredDimension; ++i)
//       {
// 	int TmpCharge = 0;
// 	int TmpCharge1;
// 	for (int j = 0; j < this->NbrSite; ++j)
// 	{
// 	  TmpCharge1 = this->FindMaximumChargeSurroundingPlaquettes(j, this->StateDescription[i]);
// 	  if (TmpCharge1 > TmpCharge)
// 	    TmpCharge = TmpCharge1;
// 	}
// 	if (TmpCharge < 7)
// 	  TmpHilbertSpaceDimension += 1;
//       }
//       cout << "filtered Hilbert space dim = " << TmpHilbertSpaceDimension << endl;
      
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  if (this->FindStateIndex(this->StateDescription[i],  this->StateHighestBit[i]) != i)
	    {
	      cout << i << " " << this->FindStateIndex(this->StateDescription[i],  this->StateHighestBit[i]) << endl;
	    }
	}


#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion(const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSitesX = fermions.NbrSitesX;
  this->NbrSitesY = fermions.NbrSitesY;  
  this->NbrSite = fermions.NbrSite;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->ListIndicesPerPlaquette = new int*[(this->NbrSite >> 1)];
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
  {
    this->ListIndicesPerPlaquette[i] = new int[6];
    for (int j = 0; j < 6; ++j)
      this->ListIndicesPerPlaquette[i][j] = fermions.ListIndicesPerPlaquette[i][j];   
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrLzValue];
  this->LargerIndicesInPlaquette = new int**[this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
    {
      this->NbrLargerIndicesInPlaquette[i][j] = fermions.NbrLargerIndicesInPlaquette[i][j];
      if (this->NbrLargerIndicesInPlaquette[i][j] != 0)
	this->LargerIndicesInPlaquette[i][j] = new int[this->NbrLargerIndicesInPlaquette[i][j]];
      else
	this->LargerIndicesInPlaquette[i][j] = 0;
      for (int k = 0; k < this->NbrLargerIndicesInPlaquette[i][j]; ++k)
	this->LargerIndicesInPlaquette[i][j][k] = fermions.LargerIndicesInPlaquette[i][j][k];
    }
  }
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::~FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion ()
{
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
    delete[] this->ListIndicesPerPlaquette[i];
  delete[] this->ListIndicesPerPlaquette;
  
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
    for (int j = 0; j < 3; ++j)
      if (this->LargerIndicesInPlaquette[i][j] != 0)
	delete[] this->LargerIndicesInPlaquette[i][j];
    
    delete[] this->NbrLargerIndicesInPlaquette[i];
    delete[] this->LargerIndicesInPlaquette[i];
  }
  delete[] this->NbrLargerIndicesInPlaquette;
  delete[] this->LargerIndicesInPlaquette;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::operator = (const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSitesX = fermions.NbrSitesX;
  this->NbrSitesY = fermions.NbrSitesY;  
  this->NbrSite = fermions.NbrSite;
  this->NbrLzValue = fermions.NbrLzValue;
  this->ListIndicesPerPlaquette = new int*[(this->NbrSite >> 1)];
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
  {
    this->ListIndicesPerPlaquette[i] = new int[6];
    for (int j = 0; j < 6; ++j)
      this->ListIndicesPerPlaquette[i][j] = fermions.ListIndicesPerPlaquette[i][j];   
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrLzValue];
  this->LargerIndicesInPlaquette = new int**[this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
    {
      this->NbrLargerIndicesInPlaquette[i][j] = fermions.NbrLargerIndicesInPlaquette[i][j];
      if (this->NbrLargerIndicesInPlaquette[i][j] != 0)
	this->LargerIndicesInPlaquette[i][j] = new int[this->NbrLargerIndicesInPlaquette[i][j]];
      else
	this->LargerIndicesInPlaquette[i][j] = 0;
      for (int k = 0; k < this->NbrLargerIndicesInPlaquette[i][j]; ++k)
	this->LargerIndicesInPlaquette[i][j][k] = fermions.LargerIndicesInPlaquette[i][j][k];
    }
  }
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::Clone()
{
  return new FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion(*this);
}

// initialize all of the arrays that will be used to implement hexagon exclusion conditions
//

void FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::InitializeHexagonArrays()
{
  this->ListIndicesPerPlaquette = new int*[this->NbrSite >> 1];
  for (int i = 0; i < this->NbrSitesX; ++i)
  {
    for (int j = 0; j < this->NbrSitesY; ++j)
    {
      int TmpIndex = this->NbrSitesY * i + j;
      this->ListIndicesPerPlaquette[TmpIndex] = new int[6];
      this->ListIndicesPerPlaquette[TmpIndex][0] = this->FindSiteIndex(0, i, j);
      this->ListIndicesPerPlaquette[TmpIndex][1] = this->FindSiteIndex(1, i, j);
      this->ListIndicesPerPlaquette[TmpIndex][2] = this->FindSiteIndex(0, i + 1, j);
      this->ListIndicesPerPlaquette[TmpIndex][3] = this->FindSiteIndex(1, i + 1, j - 1);
      this->ListIndicesPerPlaquette[TmpIndex][4] = this->FindSiteIndex(0, i + 1, j - 1);
      this->ListIndicesPerPlaquette[TmpIndex][5] = this->FindSiteIndex(1, i, j - 1);
      
    }
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrLzValue];
  this->LargerIndicesInPlaquette = new int**[this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
      this->NbrLargerIndicesInPlaquette[i][j] = 0;
    
    int j = 0;
    bool tmpFlag;
    int tmp = 0;
    while (tmp < 3)
    {
      int tmpNbrIndex = 0;
      tmpFlag = false;
      for (int k = 0; k < 6; ++k)
      {
	if (this->ListIndicesPerPlaquette[j][k] > i)
	  tmpNbrIndex += 1;
	if (this->ListIndicesPerPlaquette[j][k] == i)
	  tmpFlag = true;
      }
      ++j;
      
      if (tmpFlag)
      {
	this->NbrLargerIndicesInPlaquette[i][tmp] = tmpNbrIndex;
	if (tmpNbrIndex > 0)
	  this->LargerIndicesInPlaquette[i][tmp] = new int[this->NbrLargerIndicesInPlaquette[i][tmp]];
	else
	  this->LargerIndicesInPlaquette[i][tmp] = 0;
	++tmp;
      }
    }
    
    j = 0;
    tmp = 0;
    while (tmp < 3)
    {
      int tmpNbrIndex = 0;
      tmpFlag = false;
      for (int k = 0; k < 6; ++k)
      {
	if ((this->ListIndicesPerPlaquette[j][k] > i) && (tmpNbrIndex < this->NbrLargerIndicesInPlaquette[i][tmp]))
	{
	  this->LargerIndicesInPlaquette[i][tmp][tmpNbrIndex] = this->ListIndicesPerPlaquette[j][k];
	  tmpNbrIndex += 1;
	}
	if (this->ListIndicesPerPlaquette[j][k] == i)
	  tmpFlag = true;
      }
      ++j;
      
      if (tmpFlag)
	++tmp;
    }
  }
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site linearized index
// currentConfiguration = configuraton at the current x coordinate
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::GenerateStates(int nbrFermions, int currentSite, unsigned long currentConfiguration, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    {
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return pos;
  if (TmpCharge < 5)
  {
    Mask =  0x3ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }
  if (TmpCharge < 6)
  {
    Mask =  0x2ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
    
    Mask =  0x1ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }  
  return this->GenerateStates(nbrFermions, currentSite - 1, currentConfiguration, pos);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// nbrSpinUp = number of fermions with spin up
// currentSite = current site linearized index
// currentConfiguration = configuraton at the current x coordinate
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::GenerateStates(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration, long pos)
{
   if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    {
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return pos;
  if (TmpCharge < 5)
  {
    Mask =  0x3ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 2, nbrSpinUp - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }
  if (TmpCharge < 6)
  {
    Mask =  0x2ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 1, nbrSpinUp - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
    
    Mask =  0x1ul << (currentSite << 1);
    TmpPos = this->GenerateStates(nbrFermions - 1, nbrSpinUp, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }  
  return this->GenerateStates(nbrFermions, nbrSpinUp, currentSite - 1, currentConfiguration, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentSite = current site linearized index
// currentConfiguration = configuraton at the current x coordinate
// return value = Hilbert space dimension

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, unsigned long currentConfiguration)
{
  if (nbrFermions == 0)
    {
      return 1l;
    }
    
  if ((currentSite < 0) || (nbrFermions < 0))
    {
      return 0l;
    }

  long TmpDimension = 0l;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return TmpDimension;
  if (TmpCharge < 5)
      TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentSite - 1, currentConfiguration | (0x3ul << (currentSite << 1)));
  if (TmpCharge < 6)
  {
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, currentConfiguration | (0x2ul << (currentSite << 1)));
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, currentConfiguration | (0x1ul << (currentSite << 1)));
  }  
  TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions, currentSite - 1, currentConfiguration);
  return TmpDimension;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// nbrSpinUp = number of fermions with spin up
// currentSite = current site linearized index
// currentConfiguration = configuraton at the current x coordinate
// return value = Hilbert space dimension

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      return 1l;
    }
    
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    {
      return 0l;
    }

  long TmpDimension = 0l;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return TmpDimension;
  if (TmpCharge < 5)
      TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, nbrSpinUp - 1, currentSite - 1, currentConfiguration | (0x3ul << (currentSite << 1)));
  if (TmpCharge < 6)
  {
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, nbrSpinUp - 1, currentSite - 1, currentConfiguration | (0x2ul << (currentSite << 1)));
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, nbrSpinUp, currentSite - 1, currentConfiguration | (0x1ul << (currentSite << 1)));
  }  
  TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions, nbrSpinUp, currentSite - 1, currentConfiguration);
  return TmpDimension;
}
  
  
// compute the charge on all three hexagons surrounding one site and find the largest one
//
// siteIndex = site index
// configuraton =hilbert space configuraton
// return value = maximum charge
  
int FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FindMaximumChargeSurroundingPlaquettes (int siteIndex, unsigned long configuration)
{
  int TmpIndex;
  int TmpCharge = 0;
  int TmpChargePlaquette;
  for (int i = 0; i < 3; ++i)
  {
    TmpChargePlaquette = 0;
    if (this->NbrLargerIndicesInPlaquette[siteIndex][i] > 2)
    {
      for (int j = 0; j < this->NbrLargerIndicesInPlaquette[siteIndex][i]; ++j)
      {
	TmpIndex = this->LargerIndicesInPlaquette[siteIndex][i][j];
	TmpChargePlaquette += ((configuration >> (TmpIndex << 1)) & 0x1ul);
	TmpChargePlaquette += ((configuration >> ((TmpIndex << 1) + 1)) & 0x1ul);
      }
      
      if (TmpChargePlaquette > 4)
	TmpChargePlaquette += ((configuration >> (siteIndex << 1)) & 0x1ul);
	TmpChargePlaquette += ((configuration >> ((siteIndex << 1) + 1)) & 0x1ul);
      
      if (TmpChargePlaquette > TmpCharge)
	TmpCharge = TmpChargePlaquette;
    }
  }
  return TmpCharge;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  cout << lzmax << endl;
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
  {
    if (this->StateDescription[PosMin] == stateDescription)
      return PosMin;
    else
      return this->HilbertSpaceDimension;
  }
}  