////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of quasiholes on sphere with spin,                //
//                       pairing and all momentum sectors                     //
//                                                                            //
//                        last modification : 26/08/2016                      //
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
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairingAllMomenta.h"
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

QuasiholeOnSphereWithSpinAndPairingAllMomenta::QuasiholeOnSphereWithSpinAndPairingAllMomenta()
{
}

// basic constructor
// 
// kExclusionPrinciple = k value of the exclusion principle
// rExclusionPrinciple = r value of the exclusion principle
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// maxMomentumTransfer = maximum momentum transfer allowed in the density operators
// directory = optional path to data files
// filePrefix = prefix for all input file (should include everything related to the statistics and the geometry)
// memory = amount of memory granted for precalculations

QuasiholeOnSphereWithSpinAndPairingAllMomenta::QuasiholeOnSphereWithSpinAndPairingAllMomenta (int kValue, int rValue, int lzMax, int totalSpin, int maxMomentumTransfer, 
											      const char* directory, const char* filePrefix, unsigned long memory)
{
  this->NbrFermions = 0;
  this->TotalLz = 0;
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
  this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown = new int*[this->NbrQuasiholeEntriesSingleLayer];
  for (int i = 0; i < this->NbrQuasiholeEntriesSingleLayer; ++i)
    {
      this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[i] = 0;
    }
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
  this->LzValueDownFullSpace = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = this->GenerateStates();

  if (this->LargeHilbertSpaceDimension > 0l)
    {   
      timeval TotalStartingTime;
      gettimeofday (&(TotalStartingTime), 0);
      cout << "SingleLayerDimension = " << SingleLayerDimension << endl;
      this->SingleLayerAdAMatrices = new RealMatrix**[this->LzMax + 1];
      this->SingleLayerAnnihilationMatrices = new RealMatrix**[this->LzMax + 1];
      this->SingleLayerFullAdAMatrices = new RealMatrix***[this->LzMax + 1];
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
		{
		  this->SingleLayerAnnihilationMatrices[OperatorLzValue][i] = new RealMatrix[(2 * this->GetMaximalLzSingleLayer(i)) + 1];
		}
	      else
		this->SingleLayerAnnihilationMatrices[OperatorLzValue][i] = 0;
	    }
	  this->SingleLayerFullAdAMatrices[OperatorLzValue] = new RealMatrix**[this->LzMax + 1];
	  for (int AnnihilationOperatorLzValue = 0; AnnihilationOperatorLzValue <= this->LzMax; ++AnnihilationOperatorLzValue)
	    {
	      this->SingleLayerFullAdAMatrices[OperatorLzValue][AnnihilationOperatorLzValue] = new RealMatrix*[TmpMaxNbrFermionsPerLayer + 1];
	      for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
		{ 
		  if ((abs(AnnihilationOperatorLzValue - OperatorLzValue) <= maxMomentumTransfer) && (AnnihilationOperatorLzValue != OperatorLzValue))
		    {
		      this->SingleLayerFullAdAMatrices[OperatorLzValue][AnnihilationOperatorLzValue][i] = new RealMatrix[(2 * this->GetMaximalLzSingleLayer(i)) + 1];
		    }
		  else
		    {
		      this->SingleLayerFullAdAMatrices[OperatorLzValue][AnnihilationOperatorLzValue][i] = 0;
		    }
		}     
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
		  if (TmpLeftNbrParticles >= 0)
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
		  for (int AnnihilationOperatorLzValue = -this->LzMax; AnnihilationOperatorLzValue <= this->LzMax; AnnihilationOperatorLzValue += 2)
		    {
		      int ShiftedAnnihilationOperatorLzValue = (AnnihilationOperatorLzValue + this->LzMax) / 2;
		      if ((abs(ShiftedAnnihilationOperatorLzValue - ShiftedOperatorLzValue) <= maxMomentumTransfer) && (abs(TmpRightLz + OperatorLzValue - AnnihilationOperatorLzValue) <= MaxRightTotalLz) 
			  && (AnnihilationOperatorLzValue != OperatorLzValue))
			{
			  if (directory != 0)
			    {
			      sprintf (TmpFileName, "%s/%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d_%d.mat", directory, filePrefix, this->KValue, this->RValue, 
				       TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue, AnnihilationOperatorLzValue);
			    }
			  else
			    {
			      sprintf (TmpFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d_%d.mat", filePrefix, this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, 
				       OperatorLzValue, AnnihilationOperatorLzValue);
			    }
			  if (this->SingleLayerFullAdAMatrices[ShiftedOperatorLzValue][ShiftedAnnihilationOperatorLzValue][TmpRightNbrParticles][(TmpRightLz + MaxRightTotalLz) / 2].ReadMatrix(TmpFileName) == false)
			    {
			      cout << "error, can't read " << TmpFileName << endl;
			      return;
			    } 		
			}      
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

QuasiholeOnSphereWithSpinAndPairingAllMomenta::QuasiholeOnSphereWithSpinAndPairingAllMomenta(const QuasiholeOnSphereWithSpinAndPairingAllMomenta& fermions)
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
  this->LzValueDownFullSpace = fermions.LzValueDownFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown = fermions.FirstIndexWithNbrParticlesUpLzValueUpLzValueDown;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->SingleLayerFullAdAMatrices = fermions.SingleLayerFullAdAMatrices;
  this->SingleLayerAnnihilationMatrices = fermions.SingleLayerAnnihilationMatrices;
  this->Error = fermions.Error;
}

// destructor
//

QuasiholeOnSphereWithSpinAndPairingAllMomenta::~QuasiholeOnSphereWithSpinAndPairingAllMomenta ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LzValueDownFullSpace;
      for (int i = 0; i < this->NbrQuasiholeEntriesSingleLayer; ++i)
	{
	  if (this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown != 0)
	    {
	      delete[] this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[i];
	    }
	}
      delete[] this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown;
      int TmpMaxNbrFermionsPerLayer = this->NbrFermionsUpMax;
      if (TmpMaxNbrFermionsPerLayer < this->NbrFermionsDownMax)
	{
	  TmpMaxNbrFermionsPerLayer = this->NbrFermionsDownMax;
	}
      for (int AnnihilationOperatorLzValue = 0; AnnihilationOperatorLzValue <= this->LzMax; ++AnnihilationOperatorLzValue)
	{
	  for (int CreationOperatorLzValue = 0; CreationOperatorLzValue <= this->LzMax; ++CreationOperatorLzValue)
	    {
	      for (int i = 0; i <= TmpMaxNbrFermionsPerLayer; ++i)
		{	      
		  if (this->SingleLayerFullAdAMatrices[CreationOperatorLzValue][AnnihilationOperatorLzValue][i] != 0)
		    delete[] this->SingleLayerFullAdAMatrices[CreationOperatorLzValue][AnnihilationOperatorLzValue][i];
		}	
	      delete[] this->SingleLayerFullAdAMatrices[CreationOperatorLzValue][AnnihilationOperatorLzValue];
	    }
	}
      for (int CreationOperatorLzValue = 0; CreationOperatorLzValue <= this->LzMax; ++CreationOperatorLzValue)
	{
	  delete[] this->SingleLayerFullAdAMatrices[CreationOperatorLzValue];
	}
      delete[] this->SingleLayerFullAdAMatrices;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

QuasiholeOnSphereWithSpinAndPairingAllMomenta& QuasiholeOnSphereWithSpinAndPairingAllMomenta::operator = (const QuasiholeOnSphereWithSpinAndPairingAllMomenta& fermions)
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
  this->LzValueDownFullSpace = fermions.LzValueDownFullSpace;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->NbrFermionsDownMin = fermions.NbrFermionsDownMin;
  this->NbrFermionsDownMax = fermions.NbrFermionsDownMax;
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown = fermions.FirstIndexWithNbrParticlesUpLzValueUpLzValueDown;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
  this->MaximalLzSingleLayer = fermions.MaximalLzSingleLayer;
  this->SingleLayerLinearIndices = fermions.SingleLayerLinearIndices;
  this->SingleLayerAdAMatrices = fermions.SingleLayerAdAMatrices;
  this->SingleLayerFullAdAMatrices = fermions.SingleLayerFullAdAMatrices;
  this->SingleLayerAnnihilationMatrices = fermions.SingleLayerAnnihilationMatrices;
  this->Error = fermions.Error;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* QuasiholeOnSphereWithSpinAndPairingAllMomenta::Clone()
{
  return new QuasiholeOnSphereWithSpinAndPairingAllMomenta(*this);
}

// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AuAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  int lz = 2 * m - this->LzMax;
  
  int NbrParticlesUpRight = this->NbrFermionUpFullSpace[index];
  int NbrParticlesUpLeft = NbrParticlesUpRight - 1;
  int LzUpRight = this->LzValueUpFullSpace[index];
  int LzUpLeft = LzUpRight - lz;  
  if ((NbrParticlesUpLeft < this->NbrFermionsUpMin) || (abs(LzUpLeft) > this->GetMaximalLzSingleLayer(NbrParticlesUpLeft)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUpLeft, LzUpLeft);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUpRight, LzUpRight);
  
  int NbrParticlesDownRight = NbrParticlesUpRight - this->TotalSpin;
  int NbrParticlesDownLeft = NbrParticlesUpLeft - this->TotalSpin;
  int LzDownRight = this->LzValueDownFullSpace[index];
  int LzDownLeft = LzDownRight - lz;
  if ((NbrParticlesDownLeft < 0) || (abs(LzDownLeft) > this->GetMaximalLzSingleLayer(NbrParticlesDownLeft)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDownLeft, LzDownLeft);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDownRight, LzDownRight);
  
  
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
	  this->SingleLayerAnnihilationMatrices[m][NbrParticlesDownRight][(LzDownRight + this->MaximalLzSingleLayer[NbrParticlesDownRight]) / 2].GetMatrixElement(TmpBetaLeftDown, BetaDown, Tmp1);
	  if (fabs(Tmp1) > this->Error)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUpLeft, NbrParticlesDownLeft, LzUpLeft, LzDownLeft, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	} 
    }
  return TmpNbrCouplingElements;    
}

// apply a_u_m a_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for the left annihilation operator
// n = index for the right annihilation operator
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AuAd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  int LzUp = 2 * m - this->LzMax;
  int LzDown = 2 * n - this->LzMax;
  
  int NbrParticlesUpRight = this->NbrFermionUpFullSpace[index];
  int NbrParticlesUpLeft = NbrParticlesUpRight - 1;
  int LzUpRight = this->LzValueUpFullSpace[index];
  int LzUpLeft = LzUpRight - LzUp;  
  if ((NbrParticlesUpLeft < this->NbrFermionsUpMin) || (abs(LzUpLeft) > this->GetMaximalLzSingleLayer(NbrParticlesUpLeft)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUpLeft, LzUpLeft);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUpRight, LzUpRight);
  
  int NbrParticlesDownRight = NbrParticlesUpRight - this->TotalSpin;
  int NbrParticlesDownLeft = NbrParticlesUpLeft - this->TotalSpin;
  int LzDownRight = this->LzValueDownFullSpace[index];
  int LzDownLeft = LzDownRight - LzDown;
  if ((NbrParticlesDownLeft < 0) || (abs(LzDownLeft) > this->GetMaximalLzSingleLayer(NbrParticlesDownLeft)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDownLeft, LzDownLeft);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDownRight, LzDownRight);
  
  
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
	  this->SingleLayerAnnihilationMatrices[n][NbrParticlesDownRight][(LzDownRight + this->MaximalLzSingleLayer[NbrParticlesDownRight]) / 2].GetMatrixElement(TmpBetaLeftDown, BetaDown, Tmp1);
	  if (fabs(Tmp1) > this->Error)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUpLeft, NbrParticlesDownLeft, LzUpLeft, LzDownLeft, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	} 
    }
  return TmpNbrCouplingElements;    
}
  
// apply a^\dagger_u_m a^\dagger_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for the left annihilation operator
// n = index for the right annihilation operator
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AduAdd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  return 0;
}
  
// apply a^\dagger_u_m a_u_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AduAu (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUp = this->LzValueUpFullSpace[index];
  int LzDown = this->LzValueDownFullSpace[index];
  int TmpIndexUp = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, NbrParticlesDown, LzUp, LzDown, 0, BetaDown);  
  NumberCouplingElements = 0;
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp]; ++i)
    {
      this->SingleLayerAdAMatrices[m][NbrParticlesUp][(LzUp + this->GetMaximalLzSingleLayer(NbrParticlesUp)) / 2].GetMatrixElement(i, BetaUp, Tmp);
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

// apply a^\dagger_d_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AddAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUp = this->LzValueUpFullSpace[index];
  int LzDown = this->LzValueDownFullSpace[index];
  int TmpIndexDown = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown); 
  int NbrStatesUp = this->NbrQuasiholesPerNPerLzSingleLayer [this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp)];
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, NbrParticlesDown, LzUp, LzDown, BetaUp, 0);  
  NumberCouplingElements = 0;
   
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown]; ++i)
    {
      this->SingleLayerAdAMatrices[m][NbrParticlesDown][(LzDown + this->MaximalLzSingleLayer[NbrParticlesDown]) / 2].GetMatrixElement(i, BetaDown, Tmp);
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

// apply a^\dagger_u_m a_u_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation  operators
// n = index for annihilation operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AduAu (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  int ShiftedCreationOperatorLzValue = m;
  int ShiftedAnnihilationOperatorLzValue = n;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUpRight = this->LzValueUpFullSpace[index];
  int LzUpLeft = this->LzValueUpFullSpace[index] + 2 * (m - n);
  if (abs(LzUpLeft) > this->GetMaximalLzSingleLayer(NbrParticlesUp))
    {
      return 0;
    }
  int LzDown = this->LzValueDownFullSpace[index];
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, NbrParticlesDown, LzUpLeft, LzDown, 0, BetaDown);  
  NumberCouplingElements = 0;
  int TmpIndexUp = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUpLeft);  
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp]; ++i)
    {
      this->SingleLayerFullAdAMatrices[m][n][NbrParticlesUp][(LzUpRight + this->GetMaximalLzSingleLayer(NbrParticlesUp)) / 2].GetMatrixElement(i, BetaUp, Tmp);
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
  
// apply a^\dagger_d_m a_d_n to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation  operators
// n = index for annihilation operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::AddAd (int index, int m, int n, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUp = this->LzValueUpFullSpace[index];
  int LzDownRight = this->LzValueDownFullSpace[index];
  int LzDownLeft = LzDownRight + 2 * (m - n);
  if (abs(LzDownLeft) > this->GetMaximalLzSingleLayer(NbrParticlesDown))
    {
      return 0;
    }
  int NbrStatesUp = this->NbrQuasiholesPerNPerLzSingleLayer [this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp)];
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, NbrParticlesDown, LzUp, LzDownLeft, BetaUp, 0);  
  NumberCouplingElements = 0;
   
  int TmpIndexDown = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDownLeft); 
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown]; ++i)
    {
      this->SingleLayerFullAdAMatrices[m][n][NbrParticlesDown][(LzDownRight + this->MaximalLzSingleLayer[NbrParticlesDown]) / 2].GetMatrixElement(i, BetaDown, Tmp);
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
  
  
// apply a^\dagger_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles for the up layer input state
// lzUp = momentum of the up layer input state

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::Adu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
}
  
// apply a^\dagger_d_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles for the down layer input state
// lzDown = momentum of the down layer input state  

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::Add (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
}
  
// apply a_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles for the up layer input state
// lzUp = momentum of the up layer input state  

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::Au (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
}
  
// apply a_d_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles for the down layer input state
// lzDown = momentum of the down layer input state

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::Ad (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
}

// apply a^\dagger_u_m a_u_m to a given state defined only in the up layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesUp = number of particles in the up layer
// lzUp = momentum of the up layer eigenstate

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::AduAu (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesUp, int lzUp)
{
}
  
// apply a^\dagger_u_m a_u_m to a given state defined only in the down layer
//
// m = index for destruction operators
// inputState = reference of the state to act on
// outputState  = reference of the state where the result will be stored
// nbrParticlesDown = number of particles in the down layer
// lzDown = momentum of the down layer eigenstate

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::AddAd (int m, RealVector& inputState, RealVector& outputState, int nbrParticlesDown, int lzDown)
{
}

// generate all states corresponding to the constraints
// 
// return value = Hilbert space dimension      

long QuasiholeOnSphereWithSpinAndPairingAllMomenta::GenerateStates()
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
	  int MaxComplementaryTotalLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
	  for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	    {
	      TmpIndex = this->GetLinearIndexSingleLayer(TmpNbrParticles, TmpLz);
	      this->FirstIndexWithNbrParticlesUpLzValueUp[TmpIndex] = TmpDimension;
	      this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[TmpIndex] = new int [MaxComplementaryTotalLz + 1];
	      for (int TmpComplementaryLz = -MaxComplementaryTotalLz; TmpComplementaryLz <= MaxComplementaryTotalLz; TmpComplementaryLz += 2)
		{
		  this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[TmpIndex][(MaxComplementaryTotalLz + TmpComplementaryLz) / 2] = TmpDimension;
		  TmpComplementaryIndex = this->GetLinearIndexSingleLayer(TmpComplementaryNbrParticles, TmpComplementaryLz);
		  TmpDimensionSector = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpComplementaryIndex];
		  if (TmpDimensionSector > this->MaximalNumberCouplingElements)
		    this->MaximalNumberCouplingElements = TmpDimensionSector;
		  TmpDimension += TmpDimensionSector; 
		  for (int i = 0; i < TmpDimensionSector; ++i)
		    {
		      this->NbrFermionUpFullSpace[TmpIndex1] = TmpNbrParticles;
		      this->LzValueUpFullSpace[TmpIndex1] = TmpLz;
		      this->LzValueDownFullSpace[TmpIndex1] = TmpComplementaryLz;
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

long QuasiholeOnSphereWithSpinAndPairingAllMomenta::EvaluateHilbertSpaceDimension()
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
	  int MaxComplementaryTotalLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
	  for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	    {	
	      for (int TmpComplementaryLz = -MaxComplementaryTotalLz; TmpComplementaryLz <= MaxComplementaryTotalLz; TmpComplementaryLz += 2)
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


// find the values of beta in each layer for a given state
//
// index = state index
// alpha = reference on the value of alpha (up spin)
// beta = reference on the value of beta (down spsin)

void QuasiholeOnSphereWithSpinAndPairingAllMomenta::FindBetaIndices(int index, int& alpha, int& beta)
{
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp  - this->TotalSpin;
  int LzValueUp = this->LzValueUpFullSpace[index];
  int LzValueDown = this->LzValueDownFullSpace[index];
  int TmpMinIndex = this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[this->GetLinearIndexSingleLayer(NbrParticlesUp, LzValueUp)][(LzValueDown + this->GetMaximalLzSingleLayer(NbrParticlesDown)) /2];
  int TmpNbr = this->NbrQuasiholesPerNPerLzSingleLayer[this->GetLinearIndexSingleLayer(NbrParticlesUp, LzValueUp)];
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

// find state index from the value of each number in one layer
//
// nbrParticlesUp = number of particles with up spin
// nbrParticlesDown = number of particles with down spin
// lzValueUp = value of the angular momentum for up spins
// lzValueDown = value of the angular momentum for down spins
// alpha = index of state with up spin
// beta = index of state with down spin
// return value = corresponding index, -1 if an error occured

int QuasiholeOnSphereWithSpinAndPairingAllMomenta::FindStateIndex(int nbrParticlesUp, int nbrParticlesDown, int lzValueUp, int lzValueDown, int alpha, int beta)
{
  int SingleLayerIndex = this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp);
  int TmpStateIndex = (this->FirstIndexWithNbrParticlesUpLzValueUpLzValueDown[SingleLayerIndex][(lzValueDown + this->GetMaximalLzSingleLayer(nbrParticlesDown)) / 2] 
		       + ((beta * this->NbrQuasiholesPerNPerLzSingleLayer[SingleLayerIndex]) + alpha));  
  return TmpStateIndex;
}


