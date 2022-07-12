////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//             where only the parity of the particle number is fixed          //
//                                                                            //
//                        last modification : 10/08/2015                      //
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
#include "HilbertSpace/FermionOnLatticeRealSpaceFixedZNParity.h"
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

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeRealSpaceFixedZNParity::FermionOnLatticeRealSpaceFixedZNParity ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateLzMax = 0;  
  this->LargeHilbertSpaceDimension = 0;
  this->ParitySector = 0;
  this->ZNValue = 0;
}

// basic constructor
// 
// nbrSite = total number of sites 
// zNValue = N of the particle number modulo N
// paritySector = "parity" of the particle number (between 0 and ZNValue - 1)

FermionOnLatticeRealSpaceFixedZNParity::FermionOnLatticeRealSpaceFixedZNParity (int nbrSite, int zNValue, int paritySector)
{  
  this->NbrFermions = 0;
  this->IncNbrFermions = 0;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->ParitySector = paritySector;
  this->ZNValue = zNValue;
  this->LargeHilbertSpaceDimension = 1l << (this->NbrSite - 1);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];      
      unsigned long TmpMax = 1ul << this->NbrSite;
      for (unsigned long i = 0ul ; i < TmpMax; ++i)
	{
	  unsigned long TmpNbrFermions = 0x0ul;
	  for (int j = 0; j < this->NbrSite; ++j)
	    {
	      TmpNbrFermions += (i >> j) & 0x1ul;
	    }
	  if ((((int) TmpNbrFermions) % this->ZNValue) == this->ParitySector)
	    {
	      this->StateDescription[i >> 1] = i;
	    }
	}
      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      int CurrentLzMax = this->NbrLzValue;
      while ((CurrentLzMax > 0) && (((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul))
	--CurrentLzMax;
      this->StateLzMax[0] = CurrentLzMax;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  CurrentLzMax = this->NbrLzValue;
	  while (((this->StateDescription[i] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	    --CurrentLzMax;
	  this->StateLzMax[i] = CurrentLzMax;
 	}
      this->GenerateSignLookUpTable();
      
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
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeRealSpaceFixedZNParity::FermionOnLatticeRealSpaceFixedZNParity(const FermionOnLatticeRealSpaceFixedZNParity& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSite = fermions.NbrSite;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->ZNValue = fermions.ZNValue;
  this->ParitySector = fermions.ParitySector;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeRealSpaceFixedZNParity::~FermionOnLatticeRealSpaceFixedZNParity ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeRealSpaceFixedZNParity& FermionOnLatticeRealSpaceFixedZNParity::operator = (const FermionOnLatticeRealSpaceFixedZNParity& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
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
  this->ZNValue = fermions.ZNValue;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSite = fermions.NbrSite;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->ParitySector = fermions.ParitySector;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeRealSpaceFixedZNParity::Clone()
{
  return new FermionOnLatticeRealSpaceFixedZNParity(*this);
}

