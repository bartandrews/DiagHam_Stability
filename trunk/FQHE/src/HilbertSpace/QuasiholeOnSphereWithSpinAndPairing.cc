////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with spin and pairing             //
//      (i.e. Sz conservation but no conservation of the particle number)     //
//                                                                            //
//                        last modification : 12/04/2016                      //
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
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
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
#include <sys/time.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;


// default constructor
//

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing()
{
}

// basic constructor
// 
// kExclusionPrinciple = k value of the exclusion principle
// rExclusionPrinciple = r value of the exclusion principle
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// directory = optional path to data files
// filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
// discardPairing = if true, do not load the matrix element required to compute the pairing term

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing (int kValue, int rValue, int totalLz, int lzMax, int totalSpin, 
									  const char* directory, const char* filePrefix, bool discardPairing)
{
  this->NbrFermions = 0;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  
  this->Error = 1.0e-15;
  this->KValue = kValue;
  this->RValue = rValue;
  this->FermionFactor = 1;
  this->NbrFermionsUpMin = 0;
  this->NbrFermionsUpMax = (this->KValue * (this->LzMax + 1 + this->RValue)) / (this->RValue + (this->KValue * this->FermionFactor)); 
  this->NbrFermionsDownMin = 0;
  this->NbrFermionsDownMax = this->NbrFermionsUpMax;
  int TmpMaxNbrFermionsPerLayer = this->NbrFermionsUpMax;
  if (TmpMaxNbrFermionsPerLayer < this->NbrFermionsDownMax)
    {
      TmpMaxNbrFermionsPerLayer = this->NbrFermionsDownMax;
    }
 
  this->MaximalLzSingleLayer = new int[TmpMaxNbrFermionsPerLayer + 1];
  this->SingleLayerLinearIndices = new int*[TmpMaxNbrFermionsPerLayer + 1];
  
  
  
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
  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
  {
    this->SingleLayerLinearIndices[i] = new int[this->MaximalLzSingleLayer[i] + 1];
    for (int j = 0; j <= this->MaximalLzSingleLayer[i]; ++j)
    {
      this->SingleLayerLinearIndices[i][j] = TmpIndex;
      ++TmpIndex;
    }
  }
  
  this->NbrQuasiholeEntriesSingleLayer = TmpIndex;
  this->NbrQuasiholesPerNPerLzSingleLayer = new int [this->NbrQuasiholeEntriesSingleLayer];
  this->FirstIndexWithNbrParticlesUpLzValueUp = new int[this->NbrQuasiholeEntriesSingleLayer];
  
  TmpIndex = 0;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
    {
      for (int i = 0; i < (this->GetMaximalLzSingleLayer(TmpNbrParticles) + 1); ++i)
	{
	  this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = 0;
	  ++TmpIndex;  
	}
    }
  
  for (int i = 0; i < NbrNonZeroElements; ++i)
  {
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
      SingleLayerDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex];
  
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
      timeval TotalStartingTime;
      gettimeofday (&(TotalStartingTime), 0);
      cout << "SingleLayerDimension = " << SingleLayerDimension << endl;
      this->SingleLayerAdAMatrices = new RealMatrix**[this->LzMax + 1];
      this->SingleLayerAnnihilationMatrices = new RealMatrix**[this->LzMax + 1];
      for (int OperatorLzValue = 0; OperatorLzValue <= this->LzMax; ++OperatorLzValue)
	{
	  this->SingleLayerAdAMatrices[OperatorLzValue] = new RealMatrix*[TmpMaxNbrFermionsPerLayer + 1];
	  this->SingleLayerAnnihilationMatrices[OperatorLzValue] = new RealMatrix*[TmpMaxNbrFermionsPerLayer + 1];
	  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
	    {      
	      if (((i >= this->NbrFermionsUpMin) && (i <= this->NbrFermionsUpMax)) || ((i >= this->NbrFermionsDownMin) && (i <= this->NbrFermionsDownMax)))
		this->SingleLayerAdAMatrices[OperatorLzValue][i] = new RealMatrix[(2 * this->GetMaximalLzSingleLayer(i)) + 1];
	      else
		this->SingleLayerAdAMatrices[OperatorLzValue][i] = 0;
		
		
	      if (((i > this->NbrFermionsUpMin) && (i <= this->NbrFermionsUpMax)) || ((i > this->NbrFermionsDownMin) && (i <= this->NbrFermionsDownMax)))
		this->SingleLayerAnnihilationMatrices[OperatorLzValue][i] = new RealMatrix[(2 * this->GetMaximalLzSingleLayer(i)) + 1];
	      else
		this->SingleLayerAnnihilationMatrices[OperatorLzValue][i] = 0;
	    }
	}
      int TmpIndexRight;
      int TmpIndexLeft;
      double TmpInteraction;
      for (int TmpRightNbrParticles = this->NbrFermionsUpMin; TmpRightNbrParticles <= this->NbrFermionsUpMax; ++TmpRightNbrParticles)
	{
	  int MaxRightTotalLz = this->GetMaximalLzSingleLayer(TmpRightNbrParticles);
	  for (int TmpRightLz = -MaxRightTotalLz; TmpRightLz <= MaxRightTotalLz; TmpRightLz += 2)
	    {
	      for (int OperatorLzValue = -this->LzMax; OperatorLzValue <= this->LzMax; OperatorLzValue += 2)
		{
		  int ShiftedOperatorLzValue = (OperatorLzValue + this->LzMax) / 2;
		  int TmpLeftNbrParticles = TmpRightNbrParticles - 1;
		  if ((TmpLeftNbrParticles >= 0) && (discardPairing == false))
		    {
		      int MaxLeftTotalLz = this->GetMaximalLzSingleLayer(TmpLeftNbrParticles);
		      int TmpLeftLz = TmpRightLz - OperatorLzValue;
		      if ((TmpLeftLz >= -MaxLeftTotalLz) && (TmpLeftLz <= MaxLeftTotalLz))
			{
			  if (directory != 0)
			    {
			      sprintf (TmpFileName, "%s/%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", directory, filePrefix, this->KValue, this->RValue, 
				       TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			    }
			  else
			    { 
			      sprintf (TmpFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", filePrefix, this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			    }
			  if (this->SingleLayerAnnihilationMatrices[ShiftedOperatorLzValue][TmpRightNbrParticles][(TmpRightLz + MaxRightTotalLz) / 2].ReadMatrix(TmpFileName) == false)
			    {
			      cout << "error, can't read " << TmpFileName << endl;
			      return;
			    } 
			}
		    }
		  		  		  	  
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
      timeval TotalEndingTime;
      gettimeofday (&(TotalEndingTime), 0);
      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
      cout << "Hilbert space generated in " << Dt << "s" << endl;
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

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing(const QuasiholeOnSphereWithSpinAndPairing& fermions)
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
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->SingleLayerAnnihilationMatrices = fermions.SingleLayerAnnihilationMatrices;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->Error = fermions.Error;
}

// destructor
//

QuasiholeOnSphereWithSpinAndPairing::~QuasiholeOnSphereWithSpinAndPairing ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int TmpMaxNbrFermionsPerLayer = this->NbrFermionsUpMax;
      if (TmpMaxNbrFermionsPerLayer < this->NbrFermionsDownMax)
	{
	  TmpMaxNbrFermionsPerLayer = this->NbrFermionsDownMax;
	}
	
     for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
	{
	  if (this->SingleLayerLinearIndices[i] != 0)
	    delete[] this->SingleLayerLinearIndices[i];
	}
      for (int OperatorLzValue = 0; OperatorLzValue <= this->LzMax; ++OperatorLzValue)
	{
	  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
	    {	      
	      if ( this->SingleLayerAdAMatrices[OperatorLzValue][i] != 0)
		delete[] this->SingleLayerAdAMatrices[OperatorLzValue][i];
	    }	
	  delete[] this->SingleLayerAdAMatrices[OperatorLzValue];
	}
      delete[] this->SingleLayerAdAMatrices;
      if (this->SingleLayerAnnihilationMatrices != 0)
      {
	for (int OperatorLzValue = 0; OperatorLzValue <= this->LzMax; ++OperatorLzValue)
	{
	  for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
	    {
	      if (this->SingleLayerAnnihilationMatrices[OperatorLzValue][i] != 0)
		delete[] this->SingleLayerAnnihilationMatrices[OperatorLzValue][i];
	    }
	  delete[] this->SingleLayerAnnihilationMatrices[OperatorLzValue];
	}
      delete[] this->SingleLayerAnnihilationMatrices;
      }
      delete[] this->NbrFermionUpFullSpace;
      delete[] this->LzValueUpFullSpace;
      delete[] this->FirstIndexWithNbrParticlesUpLzValueUp;
      delete[] this->MaximalLzSingleLayer;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

QuasiholeOnSphereWithSpinAndPairing& QuasiholeOnSphereWithSpinAndPairing::operator = (const QuasiholeOnSphereWithSpinAndPairing& fermions)
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
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->SingleLayerAnnihilationMatrices = fermions.SingleLayerAnnihilationMatrices;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->Error = fermions.Error;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* QuasiholeOnSphereWithSpinAndPairing::Clone()
{
  return new QuasiholeOnSphereWithSpinAndPairing(*this);
}


  
// generate all states corresponding to the constraints
// 
// return value = Hilbert space dimension      

long QuasiholeOnSphereWithSpinAndPairing::GenerateStates()
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
      if ((TmpComplementaryNbrParticles >= 0) && (TmpComplementaryNbrParticles <= this->NbrFermionsUpMax))
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

long QuasiholeOnSphereWithSpinAndPairing::EvaluateHilbertSpaceDimension()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
    {
      int MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
      TmpComplementaryNbrParticles = TmpNbrParticles - this->TotalSpin;
      if ((TmpComplementaryNbrParticles >= 0) && (TmpComplementaryNbrParticles <= this->NbrFermionsUpMax))
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
int QuasiholeOnSphereWithSpinAndPairing::FindStateIndex(int nbrParticlesUp, int lzValueUp, int alpha, int beta)
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

void QuasiholeOnSphereWithSpinAndPairing::FindBetaIndices(int index, int& alpha, int& beta)
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

// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairing::AuAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  int lz = 2 * m - this->LzMax;
  
  int NbrParticlesUpRight = this->NbrFermionUpFullSpace[index];
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index] - 1;
  int LzUpRight = this->LzValueUpFullSpace[index];
  int LzUp = this->LzValueUpFullSpace[index] - lz;  
  if ((NbrParticlesUp < this->NbrFermionsUpMin) || (abs(LzUp) > this->GetMaximalLzSingleLayer(NbrParticlesUp)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUp + 1, LzUp + lz);
  
  int NbrParticlesDownRight = NbrParticlesUpRight - this->TotalSpin;
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzDownRight = LzUpRight - this->TotalLz;
  int LzDown = LzUp - this->TotalLz;
  if ((NbrParticlesDown < 0) || (abs(LzDown) > this->GetMaximalLzSingleLayer(NbrParticlesDown)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDown + 1, LzDown + lz);
  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  double Tmp1;
  
  int TmpNbrCouplingElements = 0;
  int TmpBetaLeftUp;
  int TmpBetaLeftDown;
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
  
  if (NumberCouplingElements == 0)
    return 0;
    
  for (int i = 0; i < NumberCouplingElements; ++i)
    {
      TmpBetaLeftUp = i / this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
      TmpBetaLeftDown = i % (this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft]);
      this->SingleLayerAnnihilationMatrices[m][NbrParticlesUpRight][(LzUpRight + this->MaximalLzSingleLayer[NbrParticlesUpRight]) / 2].GetMatrixElement(TmpBetaLeftUp, BetaUp, Tmp);
      
      if (fabs(Tmp) > this->Error)
	{
	  this->SingleLayerAnnihilationMatrices[m][NbrParticlesDownRight][(LzDownRight + this->MaximalLzSingleLayer[NbrParticlesDown + 1]) / 2].GetMatrixElement(TmpBetaLeftDown, BetaDown, Tmp1);
	  if (fabs(Tmp1) > this->Error)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUp, LzUp, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	} 
    }
  return TmpNbrCouplingElements;    
}


// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairing::AduAdd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  int lz = 2 * m - this->LzMax;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index] + 1;
  int LzUp = this->LzValueUpFullSpace[index] + lz;  
  if ((NbrParticlesUp > this->NbrFermionsUpMax) || (abs(LzUp) > this->GetMaximalLzSingleLayer(NbrParticlesUp)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUp - 1, LzUp - lz);
  
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzDown = LzUp - this->TotalLz;
  if ((NbrParticlesDown < 0) || (abs(LzDown) > this->GetMaximalLzSingleLayer(NbrParticlesDown)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDown - 1, LzDown - lz);
  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  double Tmp1;
  int TmpBetaLeftUp;
  int TmpBetaLeftDown;
  int TmpNbrCouplingElements = 0;
  
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
  
  if (NumberCouplingElements == 0)
    return 0;
  
  for (int i = 0; i < NumberCouplingElements; ++i)
    {
      TmpBetaLeftUp = i / this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
      TmpBetaLeftDown = i % (this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft]);
      this->SingleLayerAnnihilationMatrices[m][NbrParticlesUp][((LzUp + this->MaximalLzSingleLayer[NbrParticlesUp]) / 2)].GetMatrixElement(BetaUp, TmpBetaLeftUp, Tmp);
      if (fabs(Tmp) > this->Error)
	{
	  this->SingleLayerAnnihilationMatrices[m][NbrParticlesDown][(LzDown + this->MaximalLzSingleLayer[NbrParticlesDown]) / 2].GetMatrixElement(BetaDown, TmpBetaLeftDown, Tmp1);
	  if (fabs(Tmp1) > this->Error)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUp, LzUp, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	}
    }

  return TmpNbrCouplingElements;    
}


// apply a^\dagger_d_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairing::AddAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  int ShiftedOperatorLzValue = m;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUp = this->LzValueUpFullSpace[index];
  int LzDown = LzUp - this->TotalLz;
  int TmpIndexDown = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown); 
  int NbrStatesUp = this->NbrQuasiholesPerNPerLzSingleLayer [this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp)];
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, LzUp, BetaUp, 0);  
  NumberCouplingElements = 0;
   
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown]; ++i)
    {
      this->SingleLayerAdAMatrices[ShiftedOperatorLzValue][NbrParticlesDown][(LzDown + this->MaximalLzSingleLayer[NbrParticlesDown]) / 2].GetMatrixElement(i, BetaDown, Tmp);
      if (fabs(Tmp) > this->Error)
	{
	  interactionElements[NumberCouplingElements] = Tmp;
	  leftIndices[NumberCouplingElements] = TmpIndex;
	  ++NumberCouplingElements;
	}
      TmpIndex += NbrStatesUp;    
    }
  
  return NumberCouplingElements;    
}


// apply a^\dagger_u_m a_u_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state
int QuasiholeOnSphereWithSpinAndPairing::AduAu (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  int ShiftedOperatorLzValue = m;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int LzUp = this->LzValueUpFullSpace[index];
  int TmpIndexUp = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, LzUp, 0, BetaDown);  
  NumberCouplingElements = 0;
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp]; ++i)
    {
      this->SingleLayerAdAMatrices[ShiftedOperatorLzValue][NbrParticlesUp][(LzUp + this->GetMaximalLzSingleLayer(NbrParticlesUp)) / 2].GetMatrixElement(i, BetaUp, Tmp);
      if (fabs(Tmp) > this->Error)
	{
	  interactionElements[NumberCouplingElements] = Tmp;
	  leftIndices[NumberCouplingElements] = TmpIndex;
	  ++NumberCouplingElements;
	}    
      ++TmpIndex;
    }    
  return NumberCouplingElements;    
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& QuasiholeOnSphereWithSpinAndPairing::PrintState (ostream& Str, int state)
{
  int NbrFermionsUp = this->NbrFermionUpFullSpace[state];
  int NbrFermionsDown = NbrFermionsUp - this->TotalSpin;
  
  int lzValueUp = this->LzValueUpFullSpace[state];
  int lzValueDown = lzValueUp - this->TotalLz;
  
  int alpha;
  int beta;
  this->FindBetaIndices(state, alpha, beta);
   
  Str << "(N=" << NbrFermionsUp << ", lz=" << lzValueUp << ", " << alpha << ") (N=" << NbrFermionsDown << ", lz=" << lzValueDown << ", " << beta << ")";
  
  return Str;
}

// apply a^\dagger_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles for the up layer input state
// lzUp = momentum of the up layer input state

void QuasiholeOnSphereWithSpinAndPairing::Adu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
  outputState.TransposeMultiply(this->SingleLayerAnnihilationMatrices[m][nbrParticlesUp - 1][(lzUp + this->GetMaximalLzSingleLayer(nbrParticlesUp - 1)) / 2], inputState);
}
  
// apply a^\dagger_d_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles for the down layer input state
// lzDown = momentum of the down layer input state

void QuasiholeOnSphereWithSpinAndPairing::Add (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesDown, int lzDown)
{
  outputState.TransposeMultiply(this->SingleLayerAnnihilationMatrices[m][nbrParticlesDown - 1][(lzDown + this->GetMaximalLzSingleLayer(nbrParticlesDown - 1)) / 2], inputState);
}
  
// apply a_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles for the up layer input state
// lzUp = momentum of the up layer input state

void QuasiholeOnSphereWithSpinAndPairing::Au (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
  if ((this->SingleLayerAnnihilationMatrices[m][nbrParticlesUp] == 0) ||
      (this->SingleLayerAnnihilationMatrices[m][nbrParticlesUp][(lzUp + this->GetMaximalLzSingleLayer(nbrParticlesUp)) / 2].GetNbrRow() == 0))
    {
      outputState.ClearVector();
    }
  else
    {
      outputState.Multiply(this->SingleLayerAnnihilationMatrices[m][nbrParticlesUp][(lzUp + this->GetMaximalLzSingleLayer(nbrParticlesUp)) / 2], inputState);
    }
}
  
// apply a_d_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles for the down layer input state
// lzDown = momentum of the down layer input state

void QuasiholeOnSphereWithSpinAndPairing::Ad (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesDown, int lzDown)
{
  if ((this->SingleLayerAnnihilationMatrices[m][nbrParticlesDown] == 0) ||
      (this->SingleLayerAnnihilationMatrices[m][nbrParticlesDown][(lzDown + this->GetMaximalLzSingleLayer(nbrParticlesDown)) / 2].GetNbrRow() == 0))
    {
      outputState.ClearVector();
    }
  else
    {
      outputState.Multiply(this->SingleLayerAnnihilationMatrices[m][nbrParticlesDown][(lzDown + this->GetMaximalLzSingleLayer(nbrParticlesDown)) / 2], inputState);
    }
}
  
// apply a^\dagger_u_m a_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles in the up layer
// lzUp = momentum of the up layer eigenstate

void QuasiholeOnSphereWithSpinAndPairing::AduAu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{  
  outputState.Multiply(this->SingleLayerAdAMatrices[m][nbrParticlesUp][(lzUp + this->GetMaximalLzSingleLayer(nbrParticlesUp)) / 2], inputState);
}
  
// apply a^\dagger_u_m a_u_m to a given state defined only in the down layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles in the down layer
// lzDown = momentum of the down layer eigenstate

void QuasiholeOnSphereWithSpinAndPairing::AddAd (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesDown, int lzDown)
{
  outputState.Multiply(this->SingleLayerAdAMatrices[m][nbrParticlesDown][(lzDown + this->GetMaximalLzSingleLayer(nbrParticlesDown)) / 2], inputState);
}
  
// convert a given state from a given  n-body basis basis to another one
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis where state is defined
// return value = converted vector

RealVector QuasiholeOnSphereWithSpinAndPairing::ConvertToNbodyBasis(RealVector& state, QuasiholeOnSphereWithSpinAndPairing* nbodyBasis)
{
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  int NbrParticlesUp;
  int LzUp;
  int BetaUp;
  int BetaDown;  
  if (nbodyBasis->LargeHilbertSpaceDimension < this->LargeHilbertSpaceDimension)
    {
      for (int i= 0; i < nbodyBasis->HilbertSpaceDimension; ++i)
	{
	  NbrParticlesUp = nbodyBasis->NbrFermionUpFullSpace[i];
	  LzUp = nbodyBasis->LzValueUpFullSpace[i];
	  nbodyBasis->FindBetaIndices(i, BetaUp, BetaDown);  
	  int TmpIndex = this->FindStateIndex(NbrParticlesUp, LzUp, BetaUp, BetaDown);
	  if (TmpIndex < this->HilbertSpaceDimension)
	    {
	      TmpVector[TmpIndex] = state[i];
	    }
	}
    }
  else
    {
      for (int i= 0; i < this->HilbertSpaceDimension; ++i)
	{
	  NbrParticlesUp = this->NbrFermionUpFullSpace[i];
	  LzUp = this->LzValueUpFullSpace[i];
	  this->FindBetaIndices(i, BetaUp, BetaDown);  	  
	  int TmpIndex = nbodyBasis->FindStateIndex(NbrParticlesUp, LzUp, BetaUp, BetaDown);
	  if (TmpIndex < nbodyBasis->HilbertSpaceDimension)
	    {
	      TmpVector[i] = state[TmpIndex];
	    }
	}
    }
  return TmpVector;
}

// create a state from two single eigenstates
//
// eigenstateUp = reference of the eigenstate for the up layer
// nbrParticlesUp = number of particles in the up layer
// lzUp = momentum of the up layer eigenstate
// eigenstateDown = reference of the eigenstate for the down layer
// nbrParticlesDown = number of particles in the down layer
// lzDown = momentum of the down layer eigenstate
// return value = state built from the tensor product of the two single layer eigenstates

RealVector QuasiholeOnSphereWithSpinAndPairing::BuildFromTwoSingleLayerEigenstates(RealVector& eigenstateUp, int nbrParticlesUp, int lzUp,
										   RealVector& eigenstateDown, int nbrParticlesDown, int lzDown)
{
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (int i = 0; i < eigenstateUp.GetVectorDimension(); ++i)
    {
      for (int j = 0; j < eigenstateDown.GetVectorDimension(); ++j)
	{
	  int TmpIndex = this->FindStateIndex(nbrParticlesUp, lzUp, i, j);
	  if (TmpIndex < this->HilbertSpaceDimension)
	    {
	      TmpVector[TmpIndex] += eigenstateUp[i] * eigenstateDown[j];
	    }
	  
	}      
    }
  return TmpVector;
}

