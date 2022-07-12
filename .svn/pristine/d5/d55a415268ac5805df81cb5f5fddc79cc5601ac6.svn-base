////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of quasiholes on sphere with spin                 //
//               (i.e. both  Sz and the particle number conservation)         //
//                                                                            //
//                        last modification : 28/05/2016                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <cmath>
#include <bitset>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;


// default constructor
//

QuasiholeOnSphereWithSpin::QuasiholeOnSphereWithSpin()
{
}

// basic constructor
// 
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// nbrParticles = number of particles
// totalSpin = twice the total spin value
// directory = optional path to data files
// filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
// memory = amount of memory granted for precalculations

QuasiholeOnSphereWithSpin::QuasiholeOnSphereWithSpin (int kValue, int rValue, int totalLz, int lzMax, int nbrParticles, int totalSpin, 
						      const char* directory, const char* filePrefix, unsigned long memory)
{
  this->NbrFermions = nbrParticles;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  
  this->Error = 1.0e-13;
  this->KValue = kValue;
  this->RValue = rValue;
  this->FermionFactor = 1;
  this->NbrFermionsUpMin = this->NbrFermionsUp;
  this->NbrFermionsUpMax = this->NbrFermionsUp;
  this->NbrFermionsDownMin = this->NbrFermionsDown;
  this->NbrFermionsDownMax = this->NbrFermionsDown;
  this->SingleLayerAnnihilationMatrices = 0;

  int TmpMaxNbrFermionsPerLayer = this->NbrFermionsUpMax;
  if (TmpMaxNbrFermionsPerLayer < this->NbrFermionsDownMax)
    {
      TmpMaxNbrFermionsPerLayer = this->NbrFermionsDownMax;
    }
  this->MaximalLzSingleLayer = new int[TmpMaxNbrFermionsPerLayer + 1];
  this->SingleLayerLinearIndices = new int*[TmpMaxNbrFermionsPerLayer + 1];
  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
    this->SingleLayerLinearIndices[i] = 0;

  char* TmpFileName;
  if (directory != 0)
    { 
      TmpFileName = new char [512 + strlen(directory) + strlen(filePrefix)];
      sprintf (TmpFileName, "%s/%s_qh_states_k_%d_r_%d_nphi_%d.dat", directory, filePrefix, this->KValue, this->RValue, this->LzMax);
    }
  else
    {
      TmpFileName = new char [512];
      sprintf (TmpFileName, "%s_qh_states_k_%d_r_%d_nphi_%d.dat", filePrefix, this->KValue, this->RValue, this->LzMax);
    }
  MultiColumnASCIIFile DimensionQuasiholeSubspaces;
  if (DimensionQuasiholeSubspaces.Parse(TmpFileName) == false)
    {
      cout << "Error: " << TmpFileName << " cannot be found" << endl;
      this->HilbertSpaceDimension = 0;
      return;
    } 
  
  int NbrNonZeroElements = DimensionQuasiholeSubspaces.GetNbrLines();
  int* NbrFermionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (1);
  int* LzSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (2);
  int* DimensionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (3);
  
  for (int i = 0; i < NbrNonZeroElements; ++i)
    {
      if (NbrFermionsSingleLayer[i] <= TmpMaxNbrFermionsPerLayer)
	{
	  this->MaximalLzSingleLayer[NbrFermionsSingleLayer[i]] = LzSingleLayer[i];
	}
    }
  for (int i = 0; i < NbrNonZeroElements; ++i)
    {
      if ((NbrFermionsSingleLayer[i] <= TmpMaxNbrFermionsPerLayer) && (this->MaximalLzSingleLayer[NbrFermionsSingleLayer[i]] < LzSingleLayer[i])) 
	{    
	  this->MaximalLzSingleLayer[NbrFermionsSingleLayer[i]] = LzSingleLayer[i];
	}
    }
  
  int TmpIndex = 0;
  this->SingleLayerLinearIndices[this->NbrFermionsUpMin] = new int[this->MaximalLzSingleLayer[this->NbrFermionsUpMin] + 1];
  for (int j = 0; j <= this->MaximalLzSingleLayer[this->NbrFermionsUpMin]; ++j)
    {
      this->SingleLayerLinearIndices[this->NbrFermionsUpMin][j] = TmpIndex;
      ++TmpIndex;
    }
  this->SingleLayerLinearIndices[this->NbrFermionsDownMin] = new int[this->MaximalLzSingleLayer[this->NbrFermionsDownMin] + 1];
  for (int j = 0; j <= this->MaximalLzSingleLayer[this->NbrFermionsDownMin]; ++j)
    {
      this->SingleLayerLinearIndices[this->NbrFermionsDownMin][j] = TmpIndex;
      ++TmpIndex;
    }

  this->NbrQuasiholeEntriesSingleLayer = TmpIndex;
  this->NbrQuasiholesPerNPerLzSingleLayer = new int [this->NbrQuasiholeEntriesSingleLayer];
//   this->SingleLayerIndices = new int [this->NbrQuasiholeEntriesSingleLayer];
  this->FirstIndexWithNbrParticlesUpLzValueUp = new int[this->NbrQuasiholeEntriesSingleLayer];

  for (int i = 0; i < this->NbrQuasiholeEntriesSingleLayer; ++i)
    {
      this->NbrQuasiholesPerNPerLzSingleLayer[i] = 0;
//       this->SingleLayerIndices[i] = 0;
    }
  
  for (int i = 0; i < NbrNonZeroElements; ++i)
    {      
      if (((NbrFermionsSingleLayer[i] >= this->NbrFermionsUpMin) && (NbrFermionsSingleLayer[i] <= this->NbrFermionsUpMax)) || 
	  ((NbrFermionsSingleLayer[i] >= this->NbrFermionsDownMin) && (NbrFermionsSingleLayer[i] <= this->NbrFermionsDownMax)))
	{
	  TmpIndex = this->GetLinearIndexSingleLayer(NbrFermionsSingleLayer[i], LzSingleLayer[i]);
	  this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = DimensionsSingleLayer[i];
	}
    }
  delete[] NbrFermionsSingleLayer;
  delete[] LzSingleLayer;
  delete[] DimensionsSingleLayer;
  
  TmpIndex = 0;
  int SingleLayerDimension = 0 ;
  for (; TmpIndex < this->NbrQuasiholeEntriesSingleLayer; ++ TmpIndex)
    {
//       this->SingleLayerIndices[TmpIndex] = SingleLayerDimension;
      SingleLayerDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex];
    }
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension();
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->NbrFermionUpFullSpace = new int [this->LargeHilbertSpaceDimension];
  this->LzValueUpFullSpace = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = this->GenerateStates();

  if (this->LargeHilbertSpaceDimension > 0l)
    {   
      this->SingleLayerAdAMatrices = new RealMatrix**[this->LzMax + 1];
      for (int OperatorLzValue = 0; OperatorLzValue <= this->LzMax; ++OperatorLzValue)
	{
	  this->SingleLayerAdAMatrices[OperatorLzValue] = new RealMatrix*[TmpMaxNbrFermionsPerLayer + 1];
	  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
	    {	      
	      if (((i >= this->NbrFermionsUpMin) && (i <= this->NbrFermionsUpMax)) || ((i >= this->NbrFermionsDownMin) && (i <= this->NbrFermionsDownMax)))
		this->SingleLayerAdAMatrices[OperatorLzValue][i] = new RealMatrix[(2 *  this->GetMaximalLzSingleLayer(i)) + 1];
	      else
		this->SingleLayerAdAMatrices[OperatorLzValue][i] = 0;
	    }
	}
      int TmpIndexRight;
      int TmpIndexLeft;
      int TmpIndexRight1;
      int TmpIndexLeft1;
      double TmpInteraction;
      for (int TmpRightNbrParticles = this->NbrFermionsUpMin; TmpRightNbrParticles <= this->NbrFermionsUpMax; ++TmpRightNbrParticles)
	{
	  int MaxRightTotalLz = this->GetMaximalLzSingleLayer(TmpRightNbrParticles);
	  int MaxRightTotalLzDown = this->GetMaximalLzSingleLayer(TmpRightNbrParticles - this->TotalSpin);	  
	  for (int TmpRightLz = -MaxRightTotalLz; TmpRightLz <= MaxRightTotalLz; TmpRightLz += 2)
	    {
	      if (abs(TmpRightLz - this->TotalLz)  <= MaxRightTotalLzDown)
		{
		  TmpIndexRight1 = this->GetLinearIndexSingleLayer(TmpRightNbrParticles, TmpRightLz);
		  for (int OperatorLzValue = -this->LzMax; OperatorLzValue <= this->LzMax; OperatorLzValue += 2)
		    {		  
		      int ShiftedOperatorLzValue = (OperatorLzValue + this->LzMax) / 2;	  
		      if (directory != 0)
			{
			  sprintf (TmpFileName, "%s/%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", directory, filePrefix, this->KValue, this->RValue, 
				   TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			}
		      else
			{
			  sprintf (TmpFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", filePrefix, this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			}
		      if (this->SingleLayerAdAMatrices[ShiftedOperatorLzValue][TmpRightNbrParticles][(TmpRightLz + MaxRightTotalLz) / 2].ReadMatrix(TmpFileName) == false)
			{
			  cout << "error, can't read " << TmpFileName << endl;
			  return;
			} 
		    }
		}
	    }
	}
      for (int TmpRightNbrParticles = this->NbrFermionsDownMin; TmpRightNbrParticles <= this->NbrFermionsDownMax; ++TmpRightNbrParticles)
	{
	  int MaxRightTotalLz = this->GetMaximalLzSingleLayer(TmpRightNbrParticles);
	  int MaxRightTotalLzUp = this->GetMaximalLzSingleLayer(this->TotalSpin + TmpRightNbrParticles);	  
	  for (int TmpRightLz = -MaxRightTotalLz; TmpRightLz <= MaxRightTotalLz; TmpRightLz += 2)
	    {
	      if (abs(this->TotalLz + TmpRightLz)  <= MaxRightTotalLzUp)
		{
		  TmpIndexRight1 = this->GetLinearIndexSingleLayer(TmpRightNbrParticles, TmpRightLz);
		  for (int OperatorLzValue = -this->LzMax; OperatorLzValue <= this->LzMax; OperatorLzValue += 2)
		    {		  
		      int ShiftedOperatorLzValue = (OperatorLzValue + this->LzMax) / 2;	 
		      if (this->SingleLayerAdAMatrices[ShiftedOperatorLzValue][TmpRightNbrParticles][(TmpRightLz + MaxRightTotalLz) / 2].GetNbrColumn() == 0) 
			{
			  if (directory != 0)
			    {
			      sprintf (TmpFileName, "%s/%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", directory, filePrefix, this->KValue, this->RValue, 
				       TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			    }
			  else
			    {
			      sprintf (TmpFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", filePrefix, this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			    }
			  if (this->SingleLayerAdAMatrices[ShiftedOperatorLzValue][TmpRightNbrParticles][(TmpRightLz + MaxRightTotalLz) / 2].ReadMatrix(TmpFileName) == false)
			    {
			      cout << "error, can't read " << TmpFileName << endl;
			      return;
			    } 
			}
		    }
		}
	    }
	}
    }

  delete[] TmpFileName;
  
  this->Flag.Initialize();

 
#ifdef __DEBUG__
//   long UsedMemory = 0;
//   UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   if (this->NbrFermions > 0)
//     UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

QuasiholeOnSphereWithSpin::QuasiholeOnSphereWithSpin(const QuasiholeOnSphereWithSpin& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->KValue = fermions.KValue;
  this->RValue = fermions.RValue;
  this->FermionFactor = fermions.FermionFactor;
  this->NbrQuasiholeEntriesSingleLayer = fermions.NbrQuasiholeEntriesSingleLayer;
  this->NbrQuasiholesPerNPerLzSingleLayer = fermions.NbrQuasiholesPerNPerLzSingleLayer;
//   this->SingleLayerIndices = fermions.SingleLayerIndices;
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->Error = fermions.Error;
}

// destructor
//

QuasiholeOnSphereWithSpin::~QuasiholeOnSphereWithSpin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

QuasiholeOnSphereWithSpin& QuasiholeOnSphereWithSpin::operator = (const QuasiholeOnSphereWithSpin& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->KValue = fermions.KValue;
  this->RValue = fermions.RValue;
  this->FermionFactor = fermions.FermionFactor;
  this->NbrQuasiholeEntriesSingleLayer = fermions.NbrQuasiholeEntriesSingleLayer;
  this->NbrQuasiholesPerNPerLzSingleLayer = fermions.NbrQuasiholesPerNPerLzSingleLayer;
//   this->SingleLayerIndices = fermions.SingleLayerIndices;
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->Error = fermions.Error;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* QuasiholeOnSphereWithSpin::Clone()
{
  return new QuasiholeOnSphereWithSpin(*this);
}


  
// generate all states corresponding to the constraints
// 
// return value = Hilbert space dimension      

long QuasiholeOnSphereWithSpin::GenerateStates()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  int TmpDimensionSector;
  this->MaximalNumberCouplingElements = 0;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  
  int TmpIndex1 = 0;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
    {
      int MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
      TmpComplementaryNbrParticles = TmpNbrParticles - this->TotalSpin;
      if ((TmpComplementaryNbrParticles >= this->NbrFermionsDownMin) && (TmpComplementaryNbrParticles <= this->NbrFermionsDownMax))
	{
	  int MaxComplementaryLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
	  for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	    {
	      TmpComplementaryLz = TmpLz - this->TotalLz;
	      if (abs(TmpComplementaryLz) <= MaxComplementaryLz)
		{
		  TmpIndex = this->GetLinearIndexSingleLayer(TmpNbrParticles, TmpLz);
		  TmpComplementaryIndex = this->GetLinearIndexSingleLayer(TmpComplementaryNbrParticles, TmpComplementaryLz);
		  TmpDimensionSector = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[ TmpComplementaryIndex];
		  if (TmpDimensionSector > this->MaximalNumberCouplingElements)
		    this->MaximalNumberCouplingElements = TmpDimensionSector;
		  this->FirstIndexWithNbrParticlesUpLzValueUp[TmpIndex] = TmpDimension;
		  TmpDimension += TmpDimensionSector; 
		  for (int i = 0; i < TmpDimensionSector; ++i)
		    {
		      this->NbrFermionUpFullSpace[TmpIndex1] = TmpNbrParticles;
		      this->LzValueUpFullSpace[TmpIndex1] = TmpLz;
		      ++TmpIndex1;
		      
		    }
		}
	    }
	}
    }
  return TmpDimension;      
}



// evaluate Hilbert space dimension
//
// return value = Hilbert space dimension

long QuasiholeOnSphereWithSpin::EvaluateHilbertSpaceDimension()
{
  long TmpDimension = 0l;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
    {
      int MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
      TmpComplementaryNbrParticles = TmpNbrParticles - this->TotalSpin;
      if ((TmpComplementaryNbrParticles >= this->NbrFermionsDownMin) && (TmpComplementaryNbrParticles <= this->NbrFermionsDownMax))
	{
	  int MaxComplementaryLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
	  for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	    {
	      TmpComplementaryLz = TmpLz - this->TotalLz;
	      if (abs(TmpComplementaryLz) <= MaxComplementaryLz)
		{
		  TmpComplementaryIndex = this->GetLinearIndexSingleLayer(TmpComplementaryNbrParticles, TmpComplementaryLz);
		  TmpIndex = this->GetLinearIndexSingleLayer(TmpNbrParticles, TmpLz);
		  TmpDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpComplementaryIndex];  
		}
	    }
	}
    }
  return TmpDimension; 
}


// find state index from the value of each number in one layer
//
// nbrParticlesUp = number of particles with up spin
// lzValueUp = value of the angular momentum for up spins
// alpha = index of state with up spin
// beta = index of state with down spin
// return value = corresponding index, -1 if an error occured
int QuasiholeOnSphereWithSpin::FindStateIndex(int nbrParticlesUp, int lzValueUp, int alpha, int beta)
{
  int SingleLayerIndex = this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp);
  int TmpStateIndex = this->FirstIndexWithNbrParticlesUpLzValueUp[SingleLayerIndex] + beta * this->NbrQuasiholesPerNPerLzSingleLayer[SingleLayerIndex] + alpha;  
  return TmpStateIndex;
}

// find the values of beta in each layer for a given state
//
// index = state index
// alpha = reference on the value of alpha (up spin)
// beta = reference on the value of beta (down spsin)

void QuasiholeOnSphereWithSpin::FindBetaIndices(int index, int& alpha, int& beta)
{
  int nbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int lzValueUp = this->LzValueUpFullSpace[index];
  int TmpMinIndex = this->FirstIndexWithNbrParticlesUpLzValueUp[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)];
  int TmpNbr = this->NbrQuasiholesPerNPerLzSingleLayer[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)];
  if (TmpNbr != 0)
    {
      alpha = (index - TmpMinIndex) % TmpNbr;
      beta = (index - TmpMinIndex) / TmpNbr;
    }
  else
    {
      alpha = 0;
      beta = index - TmpMinIndex;
    }
}


