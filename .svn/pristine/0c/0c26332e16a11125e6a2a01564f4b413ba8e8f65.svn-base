////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of MPS matrix built from a fixed charge sector of another MPS   //
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


#ifndef FQHEMPSFIXEDQSECTORMATRIX_H
#define FQHEMPSFIXEDQSECTORMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"


class FQHEMPSFixedQSectorMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // pointer to the full MPS matrix without a fixed charge sector
  AbstractFQHEMPSMatrix* MPSMatrix;

  // number of CFT sectors
  int NbrCFTSectors;

  int BMatrixGroupSize;
  // periodicity of the charge sector
  int QPeriodicity;
  // Q sector that has to be selected
  int QSector;
  // CFT sector that has to be selected
  int CFTSector;

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
  // array that gives the index of each entry of the fixed Q B matrix within the original B matrix
  int* GlobalIndices;

  // indices to keep when performing a trace for the torus geometry
  int* TopologicalSectorIndices;
  // number of indices to keep when performing a trace for the torus geometry
  int TopologicalSectorNbrIndices;

  // number of orbitals covered by the site-dependent matrices
  int SiteDependentMatrixNbrOrbitals;
  // site dependent MPS
  SparseRealMatrix** SiteDependentMatrices;
  // signed orbital indices associated to each site-dependent matrix
  int* SiteDependentMatrixOrbitalIndices;
  //  number of site-dependent matrices for each orbital
  int* NbrSiteDependentMatrices;
  // physical indices for each site-dependent matrix
  unsigned long** SiteDependentPhysicalIndices;
  
 public:
  
  // default constructor 
  //
  FQHEMPSFixedQSectorMatrix();

  // constructor from two MPS matrices (the number of B matrices has to be identical for all of them)
  //
  // matrix = MPS matrix
  // qSector = Q sector that has to be selected (from 0 to qPeriodicity-1)
  // qPeriodicity = periodicity of the charge sector, if set to zero, guess it from the filling factor
  FQHEMPSFixedQSectorMatrix(AbstractFQHEMPSMatrix* matrix, int qSector = 0, int qPeriodicity = 0);

  // destructor
  //
  ~FQHEMPSFixedQSectorMatrix();
  
  // get the number of orbitals that associated to a set of B matrices
  //
  // return value = number oforbitals
  virtual int GetNbrOrbitals();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // create the B matrices for the block state
  //
  virtual void CreateBMatrices ();

  // get the number of CFT sectors invloved on the MPS
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

  // get a given physical index at a given orbital for the site-dependent MPS
  //
  // orbitalIndex = orbital index 
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  virtual void GetSiteDependentPhysicalIndex(int orbitalIndex, int index, unsigned long* configuration);

  // compute the site-dependent matrices
  //
  // initialOrbitalIndex = index of the first orbital
  // lastOrbitalIndex = index of the last orbital
  virtual void ComputeSiteDependentMatrices(int initialOrbitalIndex, int lastOrbitalIndex);

  // get the site-dependent matrices (real version) computed through ComputeSiteDependentMatrices
  //
  // siteDependentMatrices = reference on the site-dependent matrices
  // nbrSiteDependentMatrices = reference on the array providing the number of site-dependent matrices per orbital
  // siteDependentMatrixOrbitalIndices = reference on the array providing the orbital indices 
  // siteDependentPhysicalIndices = reference on the array providing the physical indices associated to each site-dependent matrix
  // return value = number of orbitals covered by the site-dependent matrices
  virtual int GetSiteDependentMatrices(SparseRealMatrix**& siteDependentMatrices, int*& nbrSiteDependentMatrices, int*& siteDependentMatrixOrbitalIndices, unsigned long**& siteDependentPhysicalIndices);

protected:
  
  // convert an occupation configuration to the array version
  //
  // occupationConfiguration = occupation configuration
  // configuration = array where the description of the physical index will be stored
  virtual void CoreGetPhysicalIndex(unsigned long occupationConfiguration, unsigned long* configuration);
  
};

// get the number of orbitals that associated to a set of B matrices
//
// return value = number oforbitals

inline int FQHEMPSFixedQSectorMatrix::GetNbrOrbitals()
{
  return this->BMatrixGroupSize;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

inline void FQHEMPSFixedQSectorMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  this->MPSMatrix->GetFillingFactor(numerator, denominator);
}

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSFixedQSectorMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return this->TransferMatrixLargestEigenvalueDegeneracy;;
}

// get the MPS truncation level
//
// return value = truncation level

inline int FQHEMPSFixedQSectorMatrix::GetTruncationLevel()
{
  return this->PLevel;
}

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int FQHEMPSFixedQSectorMatrix::GetNbrCFTSectors()
{
  return this->NbrCFTSectors;
}

// get the auxiliary space indices that are related to a given topological scetor
//
// topologicalSector = index of the topological sector to select
// nbrIndices = reference on the integer that will be set to the number of indices
// return value = array that contains the auxiliary space indices related to the selected topological sector

inline int* FQHEMPSFixedQSectorMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
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

inline int FQHEMPSFixedQSectorMatrix::GetMaximumOccupation()
{ 
  return this->MPSMatrix->GetMaximumOccupation();
}

// get the parent MPS matrices if the current MPS matrices have ones
//
// return value = pointer to the parent MPS matrices

inline AbstractFQHEMPSMatrix* FQHEMPSFixedQSectorMatrix::GetParentMPSMatrices()
{
  return this->MPSMatrix;
}

// get the array that gives the index of each entry of the current B matrix within its parent B matrix
//
// return value = index mapping array

inline int* FQHEMPSFixedQSectorMatrix::GetIndexMappingArray()
{
  return this->GlobalIndices;
}

#endif
