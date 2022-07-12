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
////////////////////////////////////////////////////////////////////////////////


#ifndef FQHEMPSFIXEDBONDDIMENSIONMATRIX_H
#define FQHEMPSFIXEDBONDDIMENSIONMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Matrix/RealMatrix.h"


class FQHEMPSFixedBondDimensionMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // pointer to the full MPS matrix that has to be truncated
  AbstractFQHEMPSMatrix* MPSMatrix;

  // number of CFT sectors
  int NbrCFTSectors;
  
  // |P| truncation level
  int PLevel;
  
  // left bond dimension of the truncated MPS
  int LeftBondDimensionTruncatedMPS;
  // right bond dimension of the truncated MPS
  int RightBondDimensionTruncatedMPS;

  // LeftAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the left of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  RealMatrix*** LeftAuxiliaryBasisRotation;
  // RightAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the right of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  RealMatrix*** RightAuxiliaryBasisRotation;
  
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
  // array that gives the index of each entry of the fixed Q B matrix within the original B matrix
  int* GlobalIndices;

  // indices to keep when performing a trace for the torus geometry
  int* TopologicalSectorIndices;
  // number of indices to keep when performing a trace for the torus geometry
  int TopologicalSectorNbrIndices;


protected:
    
  // compute the truncated B matrices
  //
  void ComputeTruncatedBMatrices();
  
  
 public:
  
  // default constructor 
  //
  FQHEMPSFixedBondDimensionMatrix();

  // constructor from MPS matrices
  //
  // matrix = MPS matrix
  // leftAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the left of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  // rightAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the right of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  FQHEMPSFixedBondDimensionMatrix(AbstractFQHEMPSMatrix* matrix, RealMatrix*** leftAuxiliaryBasisRotation, RealMatrix*** rightAuxiliaryBasisRotation);
  
  // constructor from a file containing the MPS matrices
  //
  //fileName = name of the file containing the MPS matrices
  // leftAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the left of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  // rightAuxiliaryBasisRotation = set of matrices converting the auxiliary space from the CFT to the Schmidt basis (to multiply on the right of the original matrix (three indices correspond to PLevel, CFTSector and QValue, respectively)
  FQHEMPSFixedBondDimensionMatrix(char* fileName, RealMatrix*** leftAuxiliaryBasisRotation, RealMatrix*** rightAuxiliaryBasisRotation);

  // destructor
  //
  ~FQHEMPSFixedBondDimensionMatrix();
  
  // get the number of orbitals that associated to a set of B matrices
  //
  // return value = number of orbitals
//   virtual int GetNbrOrbitals();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // create the B matrices
  //
  virtual void CreateBMatrices ();

  // get the number of CFT sectors involved on the MPS
  //
  // return value = number of CFT sectors
  virtual int GetNbrCFTSectors();

  // get the auxiliary space indices that are related to a given topological scetor
  //
  // topologicalSector = index of the topological sector to select
  // nbrIndices = reference on the integer that will be set to the number of indices
  // return value = array that contains the auxiliary space indices related to the selected topological sector
  virtual int* GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices);
  
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
//   virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

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

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // get the parent MPS matrices if the current MPS matrices have ones
  //
  // return value = pointer to the parent MPS matrices
  virtual AbstractFQHEMPSMatrix* GetParentMPSMatrices();

  // get the array that gives the index of each entry of the current B matrix within its parent B matrix
  //
  // return value = index mapping array
  virtual int* GetIndexMappingArray();

  // get a given physical indiex
  //
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  virtual void GetPhysicalIndex(int index, unsigned long* configuration);

};



// get the MPS truncation level
//
// return value = truncation level

inline int FQHEMPSFixedBondDimensionMatrix::GetTruncationLevel()
{
  return this->PLevel;
}

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int FQHEMPSFixedBondDimensionMatrix::GetNbrCFTSectors()
{
  return this->NbrCFTSectors;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

inline void FQHEMPSFixedBondDimensionMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  this->MPSMatrix->GetFillingFactor(numerator, denominator);
}


// get the auxiliary space indices that are related to a given topological scetor
//
// topologicalSector = index of the topological sector to select
// nbrIndices = reference on the integer that will be set to the number of indices
// return value = array that contains the auxiliary space indices related to the selected topological sector

inline int* FQHEMPSFixedBondDimensionMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
{
  nbrIndices = this->TopologicalSectorNbrIndices;
  int* TmpIndices = new int[this->TopologicalSectorNbrIndices];
  for (int i = 0; i < this->TopologicalSectorNbrIndices; ++i)
    {
      TmpIndices[i] = this->TopologicalSectorIndices[i];
    }
  return TmpIndices;
}
  
// get the maximum occupation per orbital
//
// return value = aximum occupation per orbital

inline int FQHEMPSFixedBondDimensionMatrix::GetMaximumOccupation()
{ 
  return this->MPSMatrix->GetMaximumOccupation();
}

// get the parent MPS matrices if the current MPS matrices have ones
//
// return value = pointer to the parent MPS matrices

inline AbstractFQHEMPSMatrix* FQHEMPSFixedBondDimensionMatrix::GetParentMPSMatrices()
{
  return this->MPSMatrix;
}

// get the array that gives the index of each entry of the current B matrix within its parent B matrix
//
// return value = index mapping array

inline int* FQHEMPSFixedBondDimensionMatrix::GetIndexMappingArray()
{
  return this->GlobalIndices;
}

#endif
