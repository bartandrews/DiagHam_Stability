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
//                        last modification : 11/12/2014                      //
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


#ifndef FQHEMPSSYMMETRIZEDSTATEMATRIX_H
#define FQHEMPSSYMMETRIZEDSTATEMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"


class FQHEMPSSymmetrizedStateMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // pointer to the full MPS matrix for the first layer
  AbstractFQHEMPSMatrix* MPSMatrix1;
  // pointer to the full MPS matrix for the second layer
  AbstractFQHEMPSMatrix* MPSMatrix2;

  // number of CFT sectors
  int NbrCFTSectors;

  // periodicity of the charge sector
  int QPeriodicity;
  // Q sector that has to be selected
  int QSector;
  // CFT sector that has to be selected
  int CFTSector;

  // true if the aligned sector is selected
  bool AlignedSectorFlag;

  // degeneracy of the transfer matrix largest eigenvalue
  int TransferMatrixLargestEigenvalueDegeneracy;

  // |P| truncation level
  int PLevel;

  // number of N (i.e. charge) values per  truncation level
  int** NbrNValuesPerPLevelCFTSector;
  // initial N (i.e. charge) value per  truncation level
  int** NInitialValuePerPLevelCFTSector;
  // last N (i.e. charge) value per  truncation level
  int** NLastValuePerPLevelCFTSector;

  // first linearized index for each truncation level, CFT sector, Q sector
  int*** StartingIndexPerPLevelCFTSectorQValue;
  // number of linearized indices for each truncation level, CFT sector, Q sector
  int*** NbrIndexPerPLevelCFTSectorQValue;

  int**** GlobalIndexMapper;

  // table to convert from the tensor product index to the reshuffled full auxiliary space indices 
  int* TensorProductIndexConvertion;

 public:
  
  // default constructor 
  //
  FQHEMPSSymmetrizedStateMatrix();

  // constructor from two MPS matrices
  //
  // matrix1 = MPS matrices that desribe the first state 
  // matrix2 = MPS matrices that desribe the second state 
  // antiSymmetrizeFlag = true if anti-symmetrization 
  // unalignedSectorFlag = true if the unalign sector has to be consider
  FQHEMPSSymmetrizedStateMatrix(AbstractFQHEMPSMatrix* matrix1, AbstractFQHEMPSMatrix* matrix2, bool antiSymmetrizeFlag, bool unalignedSectorFlag);

  // destructor
  //
  ~FQHEMPSSymmetrizedStateMatrix();
  
  // create the B matrices for the block state
  //
  virtual void CreateBMatrices ();

  // get the number of CFT sectors invloved on the MPS
  //
  // return value = number of CFT sectors
  virtual int GetNbrCFTSectors();

  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // get the MPS truncation level
  //
  // return value = truncation level
  virtual int GetTruncationLevel();

  // get the range for the bond index when fixing the tuncation level and the charge index
  //
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = range for the bond index with fixed tuncation level and charge index
  virtual int GetBondIndexRange(int pLevel, int qValue);

  // get the range for the bond index when fixing the tuncation level, charge and CFT sector index
  //
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // cftSector = CFT sector index of the block
  // return value = range for the bond index with fixed tuncation level and charge index
  virtual int GetBondIndexRange(int pLevel, int qValue, int cftSector);

  // get the bond index for a fixed truncation level and the charge index 
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue);

  // get the bond index for a fixed truncation level, charge and CFT sector index
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // cftSector = CFT sector index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector);

  // get the charge index range
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int& minQ, int& maxQ);

  // get the charge index range at a given truncation level and in a given CFT sector
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ);

  // compute the level and the charge index of a given matrix index
  //
  // index = matrix index
  // pLevel = reference on the level
  // qValue = reference on the charge index
  virtual void GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // print a given state of the auxiliary space
  //
  // str = reference on the output stream
  // index = index of the state
  // return value = reference on the output stream
  virtual ostream& PrintAuxiliarySpaceState(ostream& str, int index);

  // get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
  //
  // cftSector = index of the CFT sector
  // return value = Q sector shift
  virtual int GetQValueCFTSectorShift(int cftSector);

};

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

inline void FQHEMPSSymmetrizedStateMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  int TmpNumerator;
  this->MPSMatrix1->GetFillingFactor(numerator, denominator);
  this->MPSMatrix2->GetFillingFactor(TmpNumerator, denominator);
  numerator += TmpNumerator;
  if (this->AlignedSectorFlag == true) 
    {
      numerator *= 2;
      denominator *= 2;
    }
}

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSSymmetrizedStateMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return this->TransferMatrixLargestEigenvalueDegeneracy;;
}

// get the MPS truncation level
//
// return value = truncation level

inline int FQHEMPSSymmetrizedStateMatrix::GetTruncationLevel()
{
  return this->PLevel;
}

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int FQHEMPSSymmetrizedStateMatrix::GetNbrCFTSectors()
{
  return this->NbrCFTSectors;
}

#endif
