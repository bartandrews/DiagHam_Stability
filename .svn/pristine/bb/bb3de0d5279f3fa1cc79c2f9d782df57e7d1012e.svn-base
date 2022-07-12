////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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
#include "HilbertSpace/FermionOnSphereDroplet.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Endian.h"
#include "Polynomial/RationalPolynomial.h"
#include "Polynomial/LongRationalPolynomial.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorSumRationalPolynomialOperation.h"

#include <math.h>
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

FermionOnSphereDroplet::FermionOnSphereDroplet()
{
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = reference on twice the momentum total value
// lzMax = maximum Lz value reached by a fermion
// nbrFluxes1 = condition for the number of fluxes in a droplet
// maxNbrParticles1 = condition for the max number of particles in a droplet
// maxNbrHoles1 = condition for the max number of holes in a droplet
// nbrFluxes2 = secondary condition for the number of fluxes in a droplet
// maxNbrParticles2 = secondary condition for the max number of particles in a droplet
// maxNbrHoles2 = secondary condition for the max number of holes in a droplet
// memory = amount of memory granted for precalculations

FermionOnSphereDroplet::FermionOnSphereDroplet (int nbrFermions, int totalLz, int lzMax, int nbrFluxes1, int maxNbrParticles1, int maxNbrHoles1, int nbrFluxes2, int maxNbrParticles2, int maxNbrHoles2, unsigned long memory)
{  
  this->NbrFluxes1 = nbrFluxes1;
  this->MaxNbrParticles1 = maxNbrParticles1;
  this->MaxNbrHoles1 = maxNbrHoles1;

  this->NbrFluxes2 = nbrFluxes2;
  this->MaxNbrParticles2 = maxNbrParticles2;
  this->MaxNbrHoles2 = maxNbrHoles2;

  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  if (this->NbrFermions > 0)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  else
    this->LargeHilbertSpaceDimension = 1l;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  if (this->NbrFermions > 0)
    {
      this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
      if ((this->StateDescription[0l] >> this->LzMax) == 0x0ul)
	{
	  int TmpLzMax = this->LzMax;
	  for  (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      this->StateLzMax[i] = TmpLzMax;
	    }
	}
    }
  else
    {
      this->StateDescription[0] = 0x0ul; 
      this->StateLzMax[0] = 0;
    }

//***********************************************************************************
  int TruncatedNewHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
   {
     int TmpNbrParticles1 = 0;
     for (int k = 0; k < this->NbrFluxes1; k++)
      TmpNbrParticles1 += this->AdA(i, k);      

     bool Test1 = true;
     if (this->NbrFluxes1 > 0)
       if ( (TmpNbrParticles1 > this->MaxNbrParticles1) || ((this->NbrFluxes1 - TmpNbrParticles1) > this->MaxNbrHoles1) )
	 Test1 = false;

     int TmpNbrParticles2 = 0;
     for (int k = 0; k < this->NbrFluxes2; k++)
      TmpNbrParticles2 += this->AdA(i, k);      

     bool Test2 = true;
     if (this->NbrFluxes2 > 0)
       if ( (TmpNbrParticles2 > this->MaxNbrParticles2) || ((this->NbrFluxes2 - TmpNbrParticles2) > this->MaxNbrHoles2) )
	 Test2 = false;

     //Gaffnian 2/5
	if (Test1 || Test2) 	
     //Needed for 3/7 state with (2,1,2) and (6,3,6) condition	
     //if ( ((Test1 == true) && (Test2==false)) || ((Test1 == false) && (Test2==true)) || ((Test1 == true) && (Test2==true)) )
     //Generic case	
     //  if (Test1 && Test2)	
         TruncatedNewHilbertSpaceDimension++;
     //else
     //  {
     //	cout<<"Reject: "; this->PrintState(cout,i); cout<<endl;
     //  }
   }

  
//  for (int i=0;i<this->HilbertSpaceDimension;++i)
//     {
//     int TmpNbrParticles = 0;
//     for (int k = 0; k < this->NbrFluxes; k++)
//      TmpNbrParticles += this->AdA(i, k); 
//     
//     if ( (TmpNbrParticles <= this->MaxNbrParticles) && ((this->NbrFluxes - TmpNbrParticles) <= this->MaxNbrHoles) )
//        { 
//            cout<<"Keep : "; this->PrintState(cout,i); cout<<"  "<<this->StateLzMax[i]<<endl;//" i "<<i<<" "<<this->FindStateIndex(this->StateDescription[i],this->StateLzMax[i])<<endl;
//         }
//       else
//         {
//            cout<<"Reject : "; this->PrintState(cout,i); cout<<"  "<<this->StateLzMax[i]<<endl;//" i "<<i<<" "<<this->FindStateIndex(this->StateDescription[i],this->StateLzMax[i])<<endl;
//         } 
//      }

  cout << "Before truncation: dim= " << this->HilbertSpaceDimension << " ";

  unsigned long* TmpStateDescription2 = new unsigned long [TruncatedNewHilbertSpaceDimension];
  int* TmpStateLzMax2 = new int [TruncatedNewHilbertSpaceDimension];
  
  TruncatedNewHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
     {

      int TmpNbrParticles1 = 0;
      for (int k = 0; k < this->NbrFluxes1; k++)
       TmpNbrParticles1 += this->AdA(i, k);      

      bool Test1 = true;
      if (this->NbrFluxes1 > 0)
        if ( (TmpNbrParticles1 > this->MaxNbrParticles1) || ((this->NbrFluxes1 - TmpNbrParticles1) > this->MaxNbrHoles1) )
	  Test1 = false;

      int TmpNbrParticles2 = 0;
      for (int k = 0; k < this->NbrFluxes2; k++)
       TmpNbrParticles2 += this->AdA(i, k);      

      bool Test2 = true;
      if (this->NbrFluxes2 > 0)
        if ( (TmpNbrParticles2 > this->MaxNbrParticles2) || ((this->NbrFluxes2 - TmpNbrParticles2) > this->MaxNbrHoles2) )
	 Test2 = false;

     //Gaffnian 2/5
	if (Test1 || Test2) 	
     //Needed for 3/7 state with (2,1,2) and (6,3,6) condition	
     //if ( ((Test1 == true) && (Test2==false)) || ((Test1 == false) && (Test2==true)) || ((Test1 == true) && (Test2==true)) )
     //Generic case
     //  if (Test1 && Test2)
         {
           TmpStateDescription2[TruncatedNewHilbertSpaceDimension] = this->StateDescription[i];
           TmpStateLzMax2[TruncatedNewHilbertSpaceDimension] = this->StateLzMax[i];
           TruncatedNewHilbertSpaceDimension++;
         } 
     }

  this->LargeHilbertSpaceDimension = TruncatedNewHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

   cout << " After truncation: dim= " << this->HilbertSpaceDimension << endl;

  //delete[] this->StateDescription;
  //delete[] this->StateLzMax;
  this->StateDescription = TmpStateDescription2;
  this->StateLzMax = TmpStateLzMax2;
//***********************************************************************************

  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);

//  for (int i=0;i<this->HilbertSpaceDimension;++i)
//     {
//            cout<<"Keep : "; this->PrintState(cout,i); cout<<"  "<<this->StateLzMax[i]<<" i "<<i<<" "<<this->TargetSpace->FindStateIndex(this->StateDescription[i],this->StateLzMax[i])<<endl;
//     }


#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
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


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereDroplet::FermionOnSphereDroplet(const FermionOnSphereDroplet& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;

  this->NbrFluxes1 = fermions.NbrFluxes1;
  this->MaxNbrParticles1 = fermions.MaxNbrParticles1;
  this->MaxNbrHoles1 = fermions.MaxNbrHoles1;

  this->NbrFluxes2 = fermions.NbrFluxes2;
  this->MaxNbrParticles2 = fermions.MaxNbrParticles2;
  this->MaxNbrHoles2 = fermions.MaxNbrHoles2;

  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
}

// destructor
//

FermionOnSphereDroplet::~FermionOnSphereDroplet ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateLzMax != 0)
	delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      if (this->LookUpTableShift != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrLzValue; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
    }
}


// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereDroplet& FermionOnSphereDroplet::operator = (const FermionOnSphereDroplet& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;

  this->NbrFluxes1 = fermions.NbrFluxes1;
  this->MaxNbrParticles1 = fermions.MaxNbrParticles1;
  this->MaxNbrHoles1 = fermions.MaxNbrHoles1;

  this->NbrFluxes2 = fermions.NbrFluxes2;
  this->MaxNbrParticles2 = fermions.MaxNbrParticles2;
  this->MaxNbrHoles2 = fermions.MaxNbrHoles2;

  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InitializeWaveFunctionEvaluation();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereDroplet::Clone()
{
  return new FermionOnSphereDroplet(*this);
}


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereDroplet::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereDroplet*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereDroplet::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereDroplet::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereDroplet::ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[nbodyBasis.FindStateIndex(this->StateDescription[i], this->StateLzMax[i])] = state[i];
  return TmpVector;
}
