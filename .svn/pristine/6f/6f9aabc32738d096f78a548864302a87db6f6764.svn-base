////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of Moore Read state wave function on disk              //
//                                                                            //
//                        last modification : 19/09/2004                      //
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
#include "Tools/FQHEWaveFunction/MooreReadOnDiskWaveFunction.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// group maximum size in bits
#define GROUP_MAXSIZE 6


// constructor
//
// nbrParticles = number of particles
// clusterSize = number of particle per cluster

MooreReadOnDiskWaveFunction::MooreReadOnDiskWaveFunction(int nbrParticles, int clusterSize)
{
  this->NbrParticles = nbrParticles;
  this->ClusterSize = clusterSize;
  this->NbrClusters = this->NbrParticles / this->ClusterSize;
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

MooreReadOnDiskWaveFunction::MooreReadOnDiskWaveFunction(const MooreReadOnDiskWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ClusterSize = function.ClusterSize;
  this->NbrClusters = function.NbrClusters;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->Flag = function.Flag;
}

// destructor
//

MooreReadOnDiskWaveFunction::~MooreReadOnDiskWaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (unsigned long i = 0; i < this->NbrPermutations; ++i)
	delete[] this->Permutations[i];
      delete[] this->Permutations;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* MooreReadOnDiskWaveFunction::Clone ()
{
  return new MooreReadOnDiskWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex MooreReadOnDiskWaveFunction::operator ()(RealVector& x)
{
  Complex* TmpZ = new Complex[this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpZ[i].Re = x[i << 1];
      TmpZ[i].Im = x[1 + (i << 1)];
    }

  Complex Value;
  Complex Tmp;
  int ReducedNbrClusters = this->NbrClusters - 1;
  int ReducedClusterSize = (this->ClusterSize - 1) * GROUP_MAXSIZE;
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    {
      unsigned long* TmpPerm = this->Permutations[i];
      Tmp = 1.0;
      for (int j = 0; j < ReducedNbrClusters; ++j)
	for (int k = j + 1; k < this->NbrClusters; ++k)
	  {
	    Tmp *= (TmpZ[TmpPerm[j] & 0x3fl] - TmpZ[TmpPerm[k] & 0x3fl]);
	    Tmp *= (TmpZ[TmpPerm[j] & 0x3fl] - TmpZ[(TmpPerm[k] >> GROUP_MAXSIZE) & 0x3fl]);
	    for (int l = GROUP_MAXSIZE; l < ReducedClusterSize; l += GROUP_MAXSIZE)
	      {
		Tmp *= (TmpZ[(TmpPerm[j] >> l) & 0x3fl] - TmpZ[(TmpPerm[k] >> l) & 0x3fl]);
		Tmp *= (TmpZ[(TmpPerm[j] >> l) & 0x3fl] - TmpZ[(TmpPerm[k] >> (l + GROUP_MAXSIZE)) & 0x3fl]);
	      }
	    Tmp *= (TmpZ[(TmpPerm[j] >> ReducedClusterSize) & 0x3fl] - TmpZ[(TmpPerm[k] >> ReducedClusterSize) & 0x3fl]);
	    Tmp *= (TmpZ[(TmpPerm[j] >> ReducedClusterSize) & 0x3fl] - TmpZ[TmpPerm[k] & 0x3fl]);
	  }
      Value += Tmp;
    }
  Value /= (double) this->NbrPermutations;
  delete[] TmpZ;
  return Value;
}


// evaluate all permutations requested for the Moore-Read state evaluation
//

void MooreReadOnDiskWaveFunction::EvaluatePermutations()
{
  unsigned long Fact = 2;
  for (unsigned long i = 3; i <= ((unsigned long) this->NbrParticles); ++i)
    Fact *= i;
  unsigned long** Perm = new unsigned long* [Fact];

  Perm[0] = new unsigned long[this->NbrClusters];
  unsigned long* TmpPerm = Perm[0];
  int Shift = 0;   
  for (int k = 0; k < this->NbrClusters; ++k) 
    {
      TmpPerm[k] = (unsigned long) 0;
      for (int i = 0; i < this->ClusterSize; ++i)    
        TmpPerm[k] |= ((unsigned long) (i + Shift)) << (i * GROUP_MAXSIZE);
      Shift += this->ClusterSize; 
    }
  for (unsigned long i = 1; i < Fact; ++i)
    {
      Perm[i] = new unsigned long[this->NbrClusters];
      for (int k = 0; k < this->NbrClusters; ++k)
        Perm[i][k] = TmpPerm[k];
    }

  int GroupId = 0;
  int PosInGroup = 0;
  for (int Pos = 0; Pos < (this->NbrParticles - 2); ++Pos)
    {
      unsigned long Step = (unsigned long) 1;
      unsigned long Lim = (unsigned long)(this->NbrParticles - Pos - 1);
      for (unsigned long i = 2; i <= Lim; ++i) 
	Step *= i;
      unsigned long j = 0;
      unsigned long Lim2 = Fact / (Step * (Lim + 1l));
      for (unsigned long k = 0; k < Lim2; ++k)
	{
	  j += Step;
	  Lim = j + Step;
	  int GroupId2 = GroupId;
	  int PosInGroup2 = PosInGroup + 1;
	  if (PosInGroup2 == this->ClusterSize)
	    {
	      PosInGroup2 = 0;
	      ++GroupId2;
	    }
	  for (int i = Pos + 1; i < this->NbrParticles; ++i)    
	    {
	      unsigned long Mask1 = ((unsigned long) 0x3f) << (PosInGroup * GROUP_MAXSIZE);
	      unsigned long Mask2 = ((unsigned long) 0x3f) << (PosInGroup2 * GROUP_MAXSIZE);
	      unsigned long NegMask1 = ~Mask1;
	      unsigned long NegMask2 = ~Mask2;
	      int Shift = (PosInGroup2 - PosInGroup) * GROUP_MAXSIZE;   
	      if (Shift > 0)
		{
		  for (; j < Lim; ++j)
		    {
		      unsigned long Tmp1 = (Perm[j][GroupId] & Mask1) << Shift;
		      unsigned long Tmp2 = (Perm[j][GroupId2] & Mask2) >> Shift;
		      Perm[j][GroupId] &= NegMask1;
		      Perm[j][GroupId] |= Tmp2;
		      Perm[j][GroupId2] &= NegMask2;
		      Perm[j][GroupId2] |= Tmp1;
		    }
		}
	      else
		{
		  Shift *= -1;
		  for (; j < Lim; ++j)
		    {
		      unsigned long Tmp1 = (Perm[j][GroupId] & Mask1) >> Shift;
		      unsigned long Tmp2 = (Perm[j][GroupId2] & Mask2) << Shift;
		      Perm[j][GroupId] &= NegMask1;
		      Perm[j][GroupId] |= Tmp2;
		      Perm[j][GroupId2] &= NegMask2;
		      Perm[j][GroupId2] |= Tmp1;
		    }
		}
	      ++PosInGroup2;
	      if (PosInGroup2 == this->ClusterSize)
		{
		  PosInGroup2 = 0;
		  ++GroupId2;
		}
	      Lim += Step;
	    } 
	}
      ++PosInGroup;
      if (PosInGroup == this->ClusterSize)
	{
	  PosInGroup = 0;
	  ++GroupId;
	}
    }
  unsigned long Mask1 = ((unsigned long) 0x3f) << (PosInGroup * GROUP_MAXSIZE);
  unsigned long Mask2 = ((unsigned long) 0x3f) << ((PosInGroup + 1) * GROUP_MAXSIZE);
  unsigned long NegMask = ~(Mask1 | Mask2);
  unsigned long Tmp;
  this->NbrPermutations = Fact;
  for (unsigned long i = 2; i < ((unsigned long) this->NbrClusters); ++i)
    this->NbrPermutations /= i;
  this->Permutations = new unsigned long* [this->NbrPermutations];
  this->NbrPermutations = 0;
  int j;
  bool TmpFlag;
  int ReducedNbrClusters = this-> NbrClusters - 1;
  for (unsigned long i = 0; i < Fact; ++i)
    {
      j = 0;
      TmpFlag = false;
      while ((j < ReducedNbrClusters) && (TmpFlag == false))
	{
	  if ((Perm[i][j] & ((unsigned long) 0x3f)) > (Perm[i][j + 1] & ((unsigned long) 0x3f))) 
	    TmpFlag = true;
	  ++j;
	}
      if (TmpFlag == false)
	{
	  this->Permutations[this->NbrPermutations] = Perm[i];
	  ++this->NbrPermutations;
	}
      else
	delete[] Perm[i];
      ++i;
      Tmp = ((Perm[i][GroupId] & Mask1) << GROUP_MAXSIZE) | ((Perm[i][GroupId] & Mask2) >> GROUP_MAXSIZE);
      Perm[i][GroupId] &= NegMask;
      Perm[i][GroupId] |= Tmp;
      j = 0;
      TmpFlag = false;
      while ((j < ReducedNbrClusters) && (TmpFlag == false))
	{
	  if ((Perm[i][j] & ((unsigned long) 0x3f)) > (Perm[i][j + 1] & ((unsigned long) 0x3f))) 
	    TmpFlag = true;
	  ++j;
	}
      if (TmpFlag == false)
	{
	  this->Permutations[this->NbrPermutations] = Perm[i];
	  ++this->NbrPermutations;
	}
      else
	delete[]  Perm[i];
    }
  delete[] Perm;
  return;
}
