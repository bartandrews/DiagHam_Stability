////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
//                                                                            //
//                                                                            //
//                 class of MPS matrix built as a 'DMRG-type'                 //
//                     truncation of another MPS matrix                       //
//                                                                            //
//                        last modification : 19/03/2016                      //
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
// /////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSFixedBondDimensionMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "GeneralTools/ArrayTools.h"

#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSFixedBondDimensionMatrix::FQHEMPSFixedBondDimensionMatrix()
{
}

// constructor from two MPS matrices (the number of B matrices has to be identical for all of them)
//
// matrix = MPS matrix
FQHEMPSFixedBondDimensionMatrix::FQHEMPSFixedBondDimensionMatrix(AbstractFQHEMPSMatrix* matrix, RealMatrix*** leftAuxiliaryBasisRotation, RealMatrix*** rightAuxiliaryBasisRotation)
{
  this->MPSMatrix = matrix;
  this->PLevel = this->MPSMatrix->GetTruncationLevel();
  this->NbrCFTSectors = this->MPSMatrix->GetNbrCFTSectors();  
//   this->TransferMatrixLargestEigenvalueDegeneracy = MPSMatrix->TransferMatrixLargestEigenvalueDegeneracy;  
   
  
  this->TorusFlag = this->MPSMatrix->IsTorus();
  this->NbrBMatrices = this->MPSMatrix->GetNbrMatrices();
  this->LeftBondDimensionTruncatedMPS = 0; 
  this->RightBondDimensionTruncatedMPS = 0; 
  
  this->LeftAuxiliaryBasisRotation = new RealMatrix**[this->PLevel + 1];
  this->RightAuxiliaryBasisRotation = new RealMatrix**[this->PLevel + 1];
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
  {
    this->LeftAuxiliaryBasisRotation[CurrentPLevel] = new RealMatrix* [this->NbrCFTSectors];
    this->RightAuxiliaryBasisRotation[CurrentPLevel] = new RealMatrix* [this->NbrCFTSectors];
    for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
    {
      int MinQValue = 0;
      int MaxQValue = 0;
      this->MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
      
      this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector] = new RealMatrix [MaxQValue - MinQValue + 1];
      this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector] = new RealMatrix [MaxQValue - MinQValue + 1];
//       
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
      {
// 	cout << leftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue] << endl;
// 	cout << rightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue] << endl;
	this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].Copy(leftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue]);
	this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].Copy(rightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue]);
// 	
	this->LeftBondDimensionTruncatedMPS += this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].GetNbrRow();
	this->RightBondDimensionTruncatedMPS += this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].GetNbrColumn();
      }
    }
  }
  
  cout << "Bond dimension on the left = " << this->LeftBondDimensionTruncatedMPS << " , Bond dimension on the right = " << this->RightBondDimensionTruncatedMPS << endl;
  
  this->ComputeTruncatedBMatrices();
  
//   unsigned long** TmpPhysicalIndices = new unsigned long*[NbrGroupBMatrices];
//   for (int i = 0; i < NbrGroupBMatrices; ++i)
//     TmpPhysicalIndices[i] = new unsigned long [this->BMatrixGroupSize];
//   cout << "grouping " << this->BMatrixGroupSize << " B matrices (" << NbrGroupBMatrices << " matrices)" << endl;
  
//   int Step = NbrGroupBMatrices / NbrBMatricesPerOrbital;
//   int TmpOrbitalIndex = this->BMatrixGroupSize - 1;
//   unsigned long* TmpPhysicalIndex = new unsigned long[this->MPSMatrix->GetNbrOrbitals()];
//   for (int i = 0; i < NbrBMatricesPerOrbital; ++i)
//     {
//       TmpSparseBMatrices[i].Copy(this->MPSMatrix->GetMatrices()[i]);
//       this->MPSMatrix->GetPhysicalIndex(i, TmpPhysicalIndex);
//       TmpPhysicalIndices[i][TmpOrbitalIndex] = TmpPhysicalIndex[0];
//     }
//   --TmpOrbitalIndex;
//   delete[] TmpPhysicalIndex;
// 
  int GroupBMatrixDimension = 0l;
  int MinQ;
  int MaxQ;
  int* QValueCFTSectorShift  = new int [this->NbrCFTSectors];
  for (int i = 0; i < this->NbrCFTSectors; ++i)
    {
      QValueCFTSectorShift[i] = this->MPSMatrix->GetQValueCFTSectorShift(i);
    }


  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrNValuesPerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = LocalMinQ;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; ++QValue)
	    {
	      GroupBMatrixDimension += this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);       
	    }
	  QValue -= 1;
	  MaxQ = QValue;
	  if (MaxQ >= MinQ)
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = (MinQ - QValueCFTSectorShift[currentCFTSector]) ;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = (MaxQ - QValueCFTSectorShift[currentCFTSector]) ;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = this->NLastValuePerPLevelCFTSector[p][currentCFTSector] - this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] + 1;
	    }
	  else
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = 0;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = -1;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = 0;
	    }
	}
    }
  int TmpBMatrixDimension = this->LeftBondDimensionTruncatedMPS;
  
  this->GlobalIndices = new int [TmpBMatrixDimension];
  for (int i = 0; i < TmpBMatrixDimension; ++i)
    this->GlobalIndices[i] = -1;
  GroupBMatrixDimension = 0;
  this->NbrIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->StartingIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->GlobalIndexMapper = new int***[this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->StartingIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->GlobalIndexMapper[p] = new int**[this->NbrCFTSectors];      
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->GlobalIndexMapper[p][currentCFTSector] = new int*[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = LocalMinQ;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; ++QValue)
	    {
	      int LocalQSector = (QValue - MinQ + QValueCFTSectorShift[currentCFTSector]);
// 	      int MaxLocalIndex = this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);
	      int MaxLocalIndex = this->LeftAuxiliaryBasisRotation[p][currentCFTSector][QValue - MinQ].GetNbrRow();
	      this->GlobalIndexMapper[p][currentCFTSector][LocalQSector] = new int[MaxLocalIndex];      
	      this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = MaxLocalIndex;
	      this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = GroupBMatrixDimension;
	      for (int i = 0; i < MaxLocalIndex; ++i)
		{
		  this->GlobalIndices[GroupBMatrixDimension] = this->MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, p, QValue, currentCFTSector);
		  this->GlobalIndexMapper[p][currentCFTSector][LocalQSector][i] = GroupBMatrixDimension;      
		  ++GroupBMatrixDimension;
		}
	    }
	}
    }
  cout << "new B matrix dimension = " << GroupBMatrixDimension << endl;
//   SparseRealMatrix* TmpSparseGroupBMatrices2 = new SparseRealMatrix[NbrGroupBMatrices];
//   this->NbrBMatrices = 0;
//   for (int i = 0; i < NbrGroupBMatrices; ++i)
//     {
//       TmpSparseGroupBMatrices2[i] = TmpSparseGroupBMatrices[i].ExtractMatrix(GroupBMatrixDimension, GroupBMatrixDimension, this->GlobalIndices, this->GlobalIndices);
//       if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
//  	++this->NbrBMatrices;
//     }

  if (this->TorusFlag == true)
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
      for (int t = 0; t < 2; ++t)
	{
	  int TmpTopologicalSectorNbrIndices = 0;
	  int* TmpTopologicalSectorIndices = this->MPSMatrix->GetTopologicalSectorIndices(t, TmpTopologicalSectorNbrIndices);
	  if (TmpTopologicalSectorNbrIndices > 0)
	    {
	      this->TopologicalSectorNbrIndices = 0;
	      this->TopologicalSectorIndices = new int[TmpTopologicalSectorNbrIndices];
	      for (int i = 0; i < TmpTopologicalSectorNbrIndices; ++i)
		{
		  int TmpPos = SearchInArray(TmpTopologicalSectorIndices[i], this->GlobalIndices, GroupBMatrixDimension);
		  if (TmpPos >= 0)
		    {
		      this->TopologicalSectorIndices[this->TopologicalSectorNbrIndices] = TmpPos;
		      ++this->TopologicalSectorNbrIndices;
		    }
		}
	      if (this->TopologicalSectorNbrIndices == 0)
		{
		  delete[] this->TopologicalSectorIndices;
		  this->TopologicalSectorIndices = 0;
		}
	      else
		{
		  t =2;
		}
	      delete[] TmpTopologicalSectorIndices;
	    }
	}
    }
  else
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
    }
  cout  << this->NbrBMatrices << " non zero B matrices" << endl;
  
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
}


// constructor from a file containing the MPS matrices 
//
// matrix = MPS matrix
FQHEMPSFixedBondDimensionMatrix::FQHEMPSFixedBondDimensionMatrix(char* fileName, RealMatrix*** leftAuxiliaryBasisRotation, RealMatrix*** rightAuxiliaryBasisRotation)
{
  this->LoadMatrices(fileName);
  this->PLevel = this->MPSMatrix->GetTruncationLevel();
  this->NbrCFTSectors = this->MPSMatrix->GetNbrCFTSectors();  
  
  this->TorusFlag = this->MPSMatrix->IsTorus();
  this->NbrBMatrices = this->MPSMatrix->GetNbrMatrices();
  this->LeftBondDimensionTruncatedMPS = 0; 
  this->RightBondDimensionTruncatedMPS = 0; 
  
  this->LeftAuxiliaryBasisRotation = new RealMatrix**[this->PLevel + 1];
  this->RightAuxiliaryBasisRotation = new RealMatrix**[this->PLevel + 1];
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
  {
    this->LeftAuxiliaryBasisRotation[CurrentPLevel] = new RealMatrix* [this->NbrCFTSectors];
    this->RightAuxiliaryBasisRotation[CurrentPLevel] = new RealMatrix* [this->NbrCFTSectors];
    for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
    {
      int MinQValue = 0;
      int MaxQValue = 0;
      this->MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
      
      this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector] = new RealMatrix [MaxQValue - MinQValue + 1];
      this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector] = new RealMatrix [MaxQValue - MinQValue + 1];
       
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
      {
	this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].Copy(leftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue]);
	this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].Copy(rightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue]);

	this->LeftBondDimensionTruncatedMPS += this->LeftAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].GetNbrRow();
	this->RightBondDimensionTruncatedMPS += this->RightAuxiliaryBasisRotation[CurrentPLevel][CurrentCFTSector][QValue - MinQValue].GetNbrColumn();
      }
    }
  }
  
  cout << "Bond dimension on the left = " << this->LeftBondDimensionTruncatedMPS << " , Bond dimension on the right = " << this->RightBondDimensionTruncatedMPS << endl;
  
  this->ComputeTruncatedBMatrices();
  
  int GroupBMatrixDimension = 0l;
  int MinQ;
  int MaxQ;
  int* QValueCFTSectorShift  = new int [this->NbrCFTSectors];
  for (int i = 0; i < this->NbrCFTSectors; ++i)
    {
      QValueCFTSectorShift[i] = this->MPSMatrix->GetQValueCFTSectorShift(i);
    }


  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrNValuesPerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = LocalMinQ;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; ++QValue)
	    {
	      GroupBMatrixDimension += this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);       
	    }
	  QValue -= 1;
	  MaxQ = QValue;
	  if (MaxQ >= MinQ)
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = (MinQ - QValueCFTSectorShift[currentCFTSector]) ;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = (MaxQ - QValueCFTSectorShift[currentCFTSector]) ;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = this->NLastValuePerPLevelCFTSector[p][currentCFTSector] - this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] + 1;
	    }
	  else
	    {
	      this->NInitialValuePerPLevelCFTSector[p][currentCFTSector] = 0;
	      this->NLastValuePerPLevelCFTSector[p][currentCFTSector] = -1;
	      this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector] = 0;
	    }
	}
    }
  int TmpBMatrixDimension = this->LeftBondDimensionTruncatedMPS;
  
  this->GlobalIndices = new int [TmpBMatrixDimension];
  for (int i = 0; i < TmpBMatrixDimension; ++i)
    this->GlobalIndices[i] = -1;
  GroupBMatrixDimension = 0;
  this->NbrIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->StartingIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
  this->GlobalIndexMapper = new int***[this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->StartingIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
      this->GlobalIndexMapper[p] = new int**[this->NbrCFTSectors];      
      for (int currentCFTSector = 0; currentCFTSector < this->NbrCFTSectors; ++currentCFTSector)
	{
	  this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector] = new int[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  this->GlobalIndexMapper[p][currentCFTSector] = new int*[this->NbrNValuesPerPLevelCFTSector[p][currentCFTSector]];
	  int LocalMinQ;
	  int LocalMaxQ;
	  this->MPSMatrix->GetChargeIndexRange(p, currentCFTSector, LocalMinQ, LocalMaxQ);
	  MinQ = LocalMinQ;
	  int QValue = MinQ;
	  MaxQ = LocalMaxQ;
	  for (; QValue <= MaxQ; ++QValue)
	    {
	      int LocalQSector = (QValue - MinQ + QValueCFTSectorShift[currentCFTSector]);
// 	      int MaxLocalIndex = this->MPSMatrix->GetBondIndexRange(p, QValue, currentCFTSector);
	      int MaxLocalIndex = this->LeftAuxiliaryBasisRotation[p][currentCFTSector][QValue - MinQ].GetNbrRow();
	      this->GlobalIndexMapper[p][currentCFTSector][LocalQSector] = new int[MaxLocalIndex];      
	      this->NbrIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = MaxLocalIndex;
	      this->StartingIndexPerPLevelCFTSectorQValue[p][currentCFTSector][LocalQSector] = GroupBMatrixDimension;
	      for (int i = 0; i < MaxLocalIndex; ++i)
		{
		  this->GlobalIndices[GroupBMatrixDimension] = this->MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, p, QValue, currentCFTSector);
		  this->GlobalIndexMapper[p][currentCFTSector][LocalQSector][i] = GroupBMatrixDimension;      
		  ++GroupBMatrixDimension;
		}
	    }
	}
    }
  cout << "new B matrix dimension = " << GroupBMatrixDimension << endl;

  if (this->TorusFlag == true)
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
      for (int t = 0; t < 2; ++t)
	{
	  int TmpTopologicalSectorNbrIndices = 0;
	  int* TmpTopologicalSectorIndices = this->MPSMatrix->GetTopologicalSectorIndices(t, TmpTopologicalSectorNbrIndices);
	  if (TmpTopologicalSectorNbrIndices > 0)
	    {
	      this->TopologicalSectorNbrIndices = 0;
	      this->TopologicalSectorIndices = new int[TmpTopologicalSectorNbrIndices];
	      for (int i = 0; i < TmpTopologicalSectorNbrIndices; ++i)
		{
		  int TmpPos = SearchInArray(TmpTopologicalSectorIndices[i], this->GlobalIndices, GroupBMatrixDimension);
		  if (TmpPos >= 0)
		    {
		      this->TopologicalSectorIndices[this->TopologicalSectorNbrIndices] = TmpPos;
		      ++this->TopologicalSectorNbrIndices;
		    }
		}
	      if (this->TopologicalSectorNbrIndices == 0)
		{
		  delete[] this->TopologicalSectorIndices;
		  this->TopologicalSectorIndices = 0;
		}
	      else
		{
		  t =2;
		}
	      delete[] TmpTopologicalSectorIndices;
	    }
	}
    }
  else
    {
      this->TopologicalSectorIndices = 0;
      this->TopologicalSectorNbrIndices = 0;
    }
  cout  << this->NbrBMatrices << " non zero B matrices" << endl;
  
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
}

  // destructor
  //

FQHEMPSFixedBondDimensionMatrix::~FQHEMPSFixedBondDimensionMatrix()
{
  if ((this->TorusFlag == true) && (this->TopologicalSectorIndices != 0))
    {
      delete[] this->TopologicalSectorIndices;
      if (this->GlobalIndices != 0)
	delete[] this->GlobalIndices;
    }
  delete[] this->LeftAuxiliaryBasisRotation;
  delete[] this->RightAuxiliaryBasisRotation;
}

// create the B matrices for the block state
//

void FQHEMPSFixedBondDimensionMatrix::CreateBMatrices ()
{
}

// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSFixedBondDimensionMatrix::GetName()
{
  char* TmpName1 = this->MPSMatrix->GetName();
  char* TmpName2 = new char [strlen(TmpName1) + 64];
  sprintf (TmpName2, "%s_truncated_chi_%d", TmpName1, this->LeftBondDimensionTruncatedMPS);
  return TmpName2;
}


// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSFixedBondDimensionMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  cout << "warning : FQHEMPSFixedBondDimensionMatrix::GetBondIndexRange should not be used" << endl;
  return -1;  
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSFixedBondDimensionMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
      (qValue > this->NLastValuePerPLevelCFTSector[pLevel][cftSector]) || (cftSector > this->NbrCFTSectors) ||  (cftSector < 0))
    return 0;
  return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]];
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSFixedBondDimensionMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  cout << "FQHEMPSFixedBondDimensionMatrix::GetBondIndexWithFixedChargeAndPLevel should not be used" << endl;
  return -1;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSFixedBondDimensionMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{  
  return this->GlobalIndexMapper[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSFixedBondDimensionMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][0];
  for (int i = 1; i < this->NbrCFTSectors; ++i)
    {
      if (this->NInitialValuePerPLevelCFTSector[pLevel][i] < minQ)
	minQ = this->NInitialValuePerPLevelCFTSector[pLevel][i];
      if (this->NLastValuePerPLevelCFTSector[pLevel][i] > maxQ)
	maxQ = this->NLastValuePerPLevelCFTSector[pLevel][i];  
    }
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSFixedBondDimensionMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSFixedBondDimensionMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  this->MPSMatrix->GetMatrixBoundaryIndices(rowIndex, columnIndex, padding);
}

// get a given physical index
//
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void FQHEMPSFixedBondDimensionMatrix::GetPhysicalIndex(int index, unsigned long* configuration)
{
  this->MPSMatrix->GetPhysicalIndex(index, configuration);
}



// compute the truncated B matrices
//
void FQHEMPSFixedBondDimensionMatrix::ComputeTruncatedBMatrices()
{
  int NbrBMatricesPerOrbital = this->MPSMatrix->GetNbrMatrices();
  
  SparseRealMatrix* TmpSparseBMatrix = new SparseRealMatrix[NbrBMatricesPerOrbital];
  TmpSparseBMatrix = this->MPSMatrix->GetMatrices();
  this->RealBMatrices = new SparseRealMatrix[NbrBMatricesPerOrbital];
  for (int i = 0; i < NbrBMatricesPerOrbital; ++i)
  {
    this->RealBMatrices[i] = SparseRealMatrix(this->LeftBondDimensionTruncatedMPS, this->RightBondDimensionTruncatedMPS);
    for (int CurrentPLevel1 = 0; CurrentPLevel1 <= PLevel; ++CurrentPLevel1)
      {
	for (int CurrentCFTSector1 = 0; CurrentCFTSector1 < NbrCFTSectors; ++CurrentCFTSector1)
	  {
	    int MinQValue1 = 0;
	    int MaxQValue1 = 0;
	    this->MPSMatrix->GetChargeIndexRange(CurrentPLevel1, CurrentCFTSector1, MinQValue1, MaxQValue1);
	    for (int QValue1 = MinQValue1; QValue1 <= MaxQValue1; ++QValue1)
	      {
		int TmpSectorLeftBondDimension = this->LeftAuxiliaryBasisRotation[CurrentPLevel1][CurrentCFTSector1][QValue1 - MinQValue1].GetNbrRow();
		int TmpFormerLeftBondDimension = this->MPSMatrix->GetBondIndexRange(CurrentPLevel1, QValue1, CurrentCFTSector1);
		for (int CurrentPLevel2 = 0; CurrentPLevel2 <= PLevel; ++CurrentPLevel2)
		{
		  for (int CurrentCFTSector2 = 0; CurrentCFTSector2 < NbrCFTSectors; ++CurrentCFTSector2)
		    {
		      int MinQValue2 = 0;
		      int MaxQValue2 = 0;
		      this->MPSMatrix->GetChargeIndexRange(CurrentPLevel2, CurrentCFTSector2, MinQValue2, MaxQValue2);
		      for (int QValue2 = MinQValue2; QValue2 <= MaxQValue2; ++QValue2)
		      {
			int TmpSectorRightBondDimension = this->RightAuxiliaryBasisRotation[CurrentPLevel2][CurrentCFTSector2][QValue2 - MinQValue2].GetNbrColumn();
			int TmpFormerRightBondDimension = this->MPSMatrix->GetBondIndexRange(CurrentPLevel2, QValue2, CurrentCFTSector2);
			for (int TmpIndexLeft = 0; TmpIndexLeft < TmpSectorLeftBondDimension; ++TmpIndexLeft)
			 {
			   for (int TmpIndexRight = 0; TmpIndexRight < TmpSectorRightBondDimension; ++TmpIndexRight)
			   {
			      double TmpCoefficient = 0.0;
			      for (int TmpIndex = 0; TmpIndex < TmpFormerLeftBondDimension; ++TmpIndex)
			      {
				double TmpCoefficient1 = 0.0;
				for (int TmpIndex1 = 0; TmpIndex1 < TmpFormerRightBondDimension; ++TmpIndex1)
				{
				  double Tmp;
				  double Tmp1;
				  TmpSparseBMatrix[i].GetMatrixElement(TmpIndex, TmpIndex1, Tmp);
				  this->RightAuxiliaryBasisRotation[CurrentPLevel2][CurrentCFTSector2][QValue2 - MinQValue2].GetMatrixElement(TmpIndex1, TmpIndexRight, Tmp1);
				  TmpCoefficient1 += Tmp * Tmp1;
				}
				double Tmp;
				this->LeftAuxiliaryBasisRotation[CurrentPLevel1][CurrentCFTSector1][QValue1 - MinQValue1].GetMatrixElement(TmpIndexLeft, TmpIndex, Tmp);
				TmpCoefficient += Tmp * TmpCoefficient1;
			      }
			      this->RealBMatrices[i].SetMatrixElement(TmpIndexLeft, TmpIndexRight, TmpCoefficient);
			    }
			 }
		      }
		    }
		   		   
		}
	      }
	  }
      }    
  }
  delete[] TmpSparseBMatrix;
}