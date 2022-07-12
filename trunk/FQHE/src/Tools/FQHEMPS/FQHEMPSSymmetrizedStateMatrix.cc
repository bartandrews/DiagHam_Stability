////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix built from symmetrized decoupled            //
//                          copies of a model state                           //
//                                                                            //
//                        last modification : 19/03/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSSymmetrizedStateMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "HilbertSpace/BosonOnDiskShort.h"

#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSSymmetrizedStateMatrix::FQHEMPSSymmetrizedStateMatrix()
{
}

// constructor from two MPS matrices
//
// matrix1 = MPS matrices that desribe the first state 
// matrix2 = MPS matrices that desribe the second state 
// antiSymmetrizeFlag = true if anti-symmetrization 
// unalignedSectorFlag = true if the unalign sector has to be consider
 
FQHEMPSSymmetrizedStateMatrix::FQHEMPSSymmetrizedStateMatrix(AbstractFQHEMPSMatrix* matrix1, AbstractFQHEMPSMatrix* matrix2, bool antiSymmetrizeFlag, bool unalignedSectorFlag)
{
  this->MPSMatrix1 = matrix1;
  this->MPSMatrix2 = matrix2;
  this->PLevel = this->MPSMatrix1->GetTruncationLevel() + this->MPSMatrix2->GetTruncationLevel();
  this->NbrCFTSectors = 1;  
  if ((antiSymmetrizeFlag == true) && (unalignedSectorFlag == false))
    {
      this->NbrCFTSectors = 2;        
    }
  this->TransferMatrixLargestEigenvalueDegeneracy = 1;  
  this->NbrBMatrices = 0;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->AlignedSectorFlag = !unalignedSectorFlag;
  if (antiSymmetrizeFlag == false)
    {
      this->NbrBMatrices = this->MPSMatrix1->GetNbrMatrices();
      if (this->NbrBMatrices > this->MPSMatrix2->GetNbrMatrices())
	this->NbrBMatrices = this->MPSMatrix2->GetNbrMatrices();
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      cout << "symmetrizing B matrices" << endl;
      this->RealBMatrices[0] = TensorProduct(this->MPSMatrix1->GetMatrices()[0], this->MPSMatrix2->GetMatrices()[0]);
      this->RealBMatrices[0] /= M_SQRT2;
      BinomialCoefficients Binomial(this->NbrBMatrices);
      for (int j = 1; j <  this->NbrBMatrices; ++j)
	{
	  this->RealBMatrices[j] = TensorProduct(this->MPSMatrix1->GetMatrices()[j], this->MPSMatrix2->GetMatrices()[0]);
	  for (int i = 1; i <= j; ++i)
	    {
	      SparseRealMatrix TmpMatrix = TensorProduct(this->MPSMatrix1->GetMatrices()[j - i], this->MPSMatrix2->GetMatrices()[i]);
	      TmpMatrix *= Binomial.GetNumericalCoefficient(j, i);
	      this->RealBMatrices[j] = this->RealBMatrices[j] + TmpMatrix;
	    }
	  this->RealBMatrices[j] /= M_SQRT2;
	}
    }
  else
    {
      if (this->MPSMatrix1->GetNbrOrbitals() == 1)
	{
	  this->NbrBMatrices = 2;
	  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
	  cout << "symmetrizing B matrices" << endl;
	  int* TmpNbrElementPerRow = new int[this->MPSMatrix1->GetMatrices()[0].GetNbrRow()];
	  for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
	    {
	      TmpNbrElementPerRow[i] = 1;
	    }
	  SparseRealMatrix SignMatrix (this->MPSMatrix1->GetMatrices()[0].GetNbrRow(), 
				       this->MPSMatrix1->GetMatrices()[0].GetNbrColumn(), TmpNbrElementPerRow);
	  int TmpP;
	  int TmpQ;
	  for (int i = 0; i < SignMatrix.GetNbrRow(); ++i)
	    {
	      this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpP, TmpQ);
	      if ((TmpQ & 1) == 0)
		SignMatrix.SetMatrixElement(i, i, -1.0);
	      else
		SignMatrix.SetMatrixElement(i, i, 1.0);
	    }
// 	  SparseRealMatrix SignMatrix2 (this->MPSMatrix2->GetMatrices()[0].GetNbrRow(), 
// 					this->MPSMatrix2->GetMatrices()[0].GetNbrColumn(), TmpNbrElementPerRow);
// 	  for (int i = 0; i < SignMatrix2.GetNbrRow(); ++i)
// 	    {
// 	      this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(i, TmpP, TmpQ);
// 	      if ((TmpQ & 1) == 0)
// 		SignMatrix2.SetMatrixElement(i, i, 1.0);
// 	      else
// 		SignMatrix2.SetMatrixElement(i, i, -1.0);
// 	    }
	  SparseRealMatrix SignedB0Matrix = Multiply(SignMatrix, this->MPSMatrix1->GetMatrices()[0]);
//	  SparseRealMatrix SignedB1Matrix = Multiply(SignMatrix2, this->MPSMatrix2->GetMatrices()[1]);
	  this->RealBMatrices[0] = TensorProduct(this->MPSMatrix1->GetMatrices()[0], this->MPSMatrix2->GetMatrices()[0]);
	  this->RealBMatrices[1] = (TensorProduct(this->MPSMatrix1->GetMatrices()[1], this->MPSMatrix2->GetMatrices()[0])
				    + TensorProduct(SignedB0Matrix, this->MPSMatrix2->GetMatrices()[1]));
	}
      else
	{
	  this->NbrBMatrices = 1 << this->MPSMatrix1->GetNbrOrbitals();
	  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
	  cout << "symmetrizing B matrices" << endl;
	}
    }
  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NbrIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  this->GlobalIndexMapper = new int*** [this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->NbrNValuesPerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[p] = new int [this->NbrCFTSectors];
      this->NbrIndexPerPLevelCFTSectorQValue[p] = new int* [this->NbrCFTSectors];
      this->GlobalIndexMapper[p] = new int** [this->NbrCFTSectors];
      for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	{
	  int MinQ = 1 << 30;
	  int MaxQ = -(1 << 30);
	  int TmpQ1 = 0;
	  int TmpQ2 = 0;
	  int TmpPLevel1 = 0;
	  int TmpPLevel2 = 0;
	  if (this->AlignedSectorFlag == true)
	    {
	      for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
		{
		  this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpPLevel1, TmpQ1);
		  for (int j = 0; j < this->MPSMatrix2->GetMatrices()[0].GetNbrRow(); ++j)
		    {
		      this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(j, TmpPLevel2, TmpQ2);
		      if (((TmpPLevel1 + TmpPLevel2 + (((TmpQ1 - TmpQ2) * (TmpQ1 - TmpQ2)) / 12)) == p) 
			  && (((abs(TmpQ1 - TmpQ2)) % 6) == (CurrentCFTSector * 3)))
			{
			  if ((TmpQ1 + TmpQ2) < MinQ)
			    {
			      MinQ = TmpQ1 + TmpQ2;
			    }
			  if ((TmpQ1 + TmpQ2) > MaxQ)
			    {
			      MaxQ = TmpQ1 + TmpQ2;
			    }
			}
		    }
		}	      
	    }
	  else
	    {
	      for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
		{
		  this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpPLevel1, TmpQ1);
		  for (int j = 0; j < this->MPSMatrix2->GetMatrices()[0].GetNbrRow(); ++j)
		    {
		      this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(j, TmpPLevel2, TmpQ2);
		      if (((TmpPLevel1 + TmpPLevel2 + (((TmpQ1 - TmpQ2) * (TmpQ1 - TmpQ2)) / 12)) == p) 
			  && (((TmpQ1 > TmpQ2) && (((TmpQ1 - TmpQ2) % 3) == 1)) 
			      || ((TmpQ1 < TmpQ2) && (((TmpQ2 - TmpQ1) % 3) == 2))))
			{
			  if ((TmpQ1 + TmpQ2) < MinQ)
			    {
			      MinQ = TmpQ1 + TmpQ2;
			    }
			  if ((TmpQ1 + TmpQ2) > MaxQ)
			    {
			      MaxQ = TmpQ1 + TmpQ2;
			    }
			}
		    }
		}
	    }
	  if (MinQ == (1 << 30))
	    {
	      this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector] = 0;
	      this->NLastValuePerPLevelCFTSector[p][CurrentCFTSector] = -1;
	      this->NbrNValuesPerPLevelCFTSector[p][CurrentCFTSector] = 0;
	      this->NbrIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector] = 0;
	      this->GlobalIndexMapper[p][CurrentCFTSector] = 0;
	    }
	  else
	    {
	      this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector] = MinQ;
	      this->NLastValuePerPLevelCFTSector[p][CurrentCFTSector] = MaxQ;
	      this->NbrNValuesPerPLevelCFTSector[p][CurrentCFTSector] = this->NLastValuePerPLevelCFTSector[p][CurrentCFTSector] - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector] + 1;
	      this->NbrIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector] = new int[MaxQ - MinQ + 1];
	      this->GlobalIndexMapper[p][CurrentCFTSector] = new int*[MaxQ - MinQ + 1];
	    }
	}
    }
  int TmpFullMatrixSize = this->MPSMatrix1->GetMatrices()[0].GetNbrRow() *  this->MPSMatrix2->GetMatrices()[0].GetNbrRow();
  this->TensorProductIndexConvertion = new int [TmpFullMatrixSize];
  for (int i = 0; i < TmpFullMatrixSize; ++i)
    this->TensorProductIndexConvertion[i] = -1;
  int TmpIndex = 0;
  this->StartingIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      this->StartingIndexPerPLevelCFTSectorQValue[p] = new int* [this->NbrCFTSectors];
      for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	{
	  this->StartingIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector] = new int [this->NbrNValuesPerPLevelCFTSector[p][CurrentCFTSector]];
	  int MinQ = this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector];
	  int MaxQ = this->NLastValuePerPLevelCFTSector[p][CurrentCFTSector];
	  int TmpQ1 = 0;
	  int TmpQ2 = 0;
	  int TmpPLevel1 = 0;
	  int TmpPLevel2 = 0;
	  for (; MinQ <= MaxQ; ++MinQ)
	    {
	      this->StartingIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]] = TmpIndex;
	      int InitialTmpIndex = TmpIndex;
	      if (this->AlignedSectorFlag == true)
		{
		  for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
		    {
		      this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpPLevel1, TmpQ1);
		      for (int j = 0; j < this->MPSMatrix2->GetMatrices()[0].GetNbrRow(); ++j)
			{
			  this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(j, TmpPLevel2, TmpQ2);
			  if (((TmpQ1 + TmpQ2) == MinQ) 
			      && ((TmpPLevel1 + TmpPLevel2 + (((TmpQ1 - TmpQ2) * (TmpQ1 - TmpQ2)) / 12)) == p) 
			      && ((abs(TmpQ1 - TmpQ2)) % 6) == (CurrentCFTSector * 3))
			    {
			      this->TensorProductIndexConvertion[(j * this->MPSMatrix1->GetMatrices()[0].GetNbrRow()) + i] = TmpIndex;
			      ++TmpIndex;
			    }
			}
		    }
		}
	      else
		{
		  for (int i = 0; i < this->MPSMatrix1->GetMatrices()[0].GetNbrRow(); ++i)
		    {
		      this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(i, TmpPLevel1, TmpQ1);
		      for (int j = 0; j < this->MPSMatrix2->GetMatrices()[0].GetNbrRow(); ++j)
			{
			  this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(j, TmpPLevel2, TmpQ2);
			  if (((TmpQ1 + TmpQ2) == MinQ) 
			      && ((TmpPLevel1 + TmpPLevel2 + (((TmpQ1 - TmpQ2) * (TmpQ1 - TmpQ2)) / 12)) == p) 
			      && (((TmpQ1 > TmpQ2) && (((TmpQ1 - TmpQ2) % 3) == 1)) 
				  || ((TmpQ1 < TmpQ2) && (((TmpQ2 - TmpQ1) % 3) == 2))))
			    {
			      this->TensorProductIndexConvertion[(j * this->MPSMatrix1->GetMatrices()[0].GetNbrRow()) + i] = TmpIndex;
			      ++TmpIndex;
			    }
			}
		    }
		}
	      this->NbrIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]] = TmpIndex - InitialTmpIndex;
	      if (TmpIndex > InitialTmpIndex)
		{
		  this->GlobalIndexMapper[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]] = new int[this->NbrIndexPerPLevelCFTSectorQValue[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]]];
		  for (int i = InitialTmpIndex; i < TmpIndex; ++i)
		    {
		      this->GlobalIndexMapper[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]][i - InitialTmpIndex] = i;
		    }
		}
	      else
		{ 
		  this->GlobalIndexMapper[p][CurrentCFTSector][MinQ - this->NInitialValuePerPLevelCFTSector[p][CurrentCFTSector]] = 0;
		}
	    }
	}
    }
  cout << "actual B matrix size " << TmpIndex << endl;
  int* TmpKeptIndices2 = new int [TmpIndex];
  for (int i = 0; i < TmpFullMatrixSize; ++i)
    {
      if (this->TensorProductIndexConvertion[i] >= 0)
	TmpKeptIndices2[this->TensorProductIndexConvertion[i]] = i;
    }
  for (int j = 0; j <  this->NbrBMatrices; ++j)
    {
      this->RealBMatrices[j] =  this->RealBMatrices[j].ExtractMatrix(TmpIndex, TmpIndex, TmpKeptIndices2, TmpKeptIndices2);
    }
  delete[] TmpKeptIndices2;

  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
}

// create the B matrices for the block state
//

void FQHEMPSSymmetrizedStateMatrix::CreateBMatrices ()
{
}

// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSSymmetrizedStateMatrix::GetName()
{
  char* TmpName1 = this->MPSMatrix1->GetName();
  char* TmpName2 = this->MPSMatrix2->GetName();
  char* TmpName3 = new char [strlen(TmpName1) + strlen(TmpName2) + 256];
  if (this->AlignedSectorFlag == true)
    {
      sprintf (TmpName3, "symmetrized_%s_%s", TmpName1, TmpName2);
    }
  else
    {
      sprintf (TmpName3, "symmetrized_unaligned_%s_%s", TmpName1, TmpName2);
    }
  return TmpName3;
}


// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  cout << "warning : FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange should not be used" << endl;
  return -1;  
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
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

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  cout << "FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargeAndPLevel should not be used" << endl;
  return -1;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSSymmetrizedStateMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{  
  return this->GlobalIndexMapper[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][localIndex];
}

// get the charge index range
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
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

void FQHEMPSSymmetrizedStateMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// compute the level and the charge index of a given matrix index
//
// index = matrix index
// pLevel = reference on the level
// qValue = reference on the charge index

void FQHEMPSSymmetrizedStateMatrix::GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue)
{
  pLevel = 0;
  for (int p = 1; p <= this->PLevel; ++p)
    {
      if (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][0][0] <= index)
	{
	  pLevel = p;
	}
    }
  
  qValue = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  for (int MinQ = 1; MinQ < this->NbrNValuesPerPLevelCFTSector[pLevel][0]; ++MinQ)
    {
      if (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][0][MinQ] <= index)
	{
	  qValue = this->NInitialValuePerPLevelCFTSector[pLevel][0] + MinQ;
	}
    }
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSSymmetrizedStateMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int TmpRowIndex1;
  int TmpColumnIndex1;
  this->MPSMatrix1->GetMatrixBoundaryIndices(TmpRowIndex1, TmpColumnIndex1, padding);
  int TmpRowIndex2;
  int TmpColumnIndex2;
  this->MPSMatrix2->GetMatrixBoundaryIndices(TmpRowIndex2, TmpColumnIndex2, padding);
  ++TmpRowIndex1;
  ++TmpColumnIndex1;
  if (this->AlignedSectorFlag == true)
    {
      ++TmpRowIndex2;
      ++TmpColumnIndex2;
    }
  rowIndex = this->TensorProductIndexConvertion[TmpRowIndex1 + this->MPSMatrix1->GetMatrices()[0].GetNbrRow() * TmpRowIndex2];
  columnIndex = this->TensorProductIndexConvertion[TmpColumnIndex1 + this->MPSMatrix1->GetMatrices()[0].GetNbrColumn() * TmpColumnIndex2];
}

// print a given state of the auxiliary space
//
// str = reference on the output stream
// index = index of the state
// return value = reference on the output stream

ostream& FQHEMPSSymmetrizedStateMatrix::PrintAuxiliarySpaceState(ostream& str, int index)
{
  int TmpPLevel1;
  int TmpQ1;
  this->MPSMatrix1->GetChargeAndPLevelFromMatrixIndex(index % this->MPSMatrix1->GetMatrices()[0].GetNbrRow(), TmpPLevel1, TmpQ1);
  int TmpPLevel2;
  int TmpQ2;
  this->MPSMatrix2->GetChargeAndPLevelFromMatrixIndex(index / this->MPSMatrix1->GetMatrices()[0].GetNbrRow(), TmpPLevel2, TmpQ2);
  str << "|" << (index  % this->MPSMatrix1->GetMatrices()[0].GetNbrRow()) << ": Q1=" << TmpQ1 << " P1=" << TmpPLevel1 
      << "> x |" << (index / this->MPSMatrix1->GetMatrices()[0].GetNbrRow()) << ": Q2=" << TmpQ2 << " P2=" << TmpPLevel2 << ">";
//   int TmpPLevel;
//   int TmpQ;
//   this->GetChargeAndPLevelFromMatrixIndex(index, TmpPLevel, TmpQ);
//   str << "|" << index << ": Q=" << TmpQ << " P=" << TmpPLevel << ">";
  return str;
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSSymmetrizedStateMatrix::GetQValueCFTSectorShift(int cftSector)
{
  if (this->AlignedSectorFlag == true) 
    {
      return (cftSector * 3);
    }
  else
    {
      return 0;
    }
}

