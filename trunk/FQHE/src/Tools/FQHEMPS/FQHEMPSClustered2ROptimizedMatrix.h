////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                                                                            //
//                        last modification : 27/09/2013                      //
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


#ifndef FQHEMPSCLUSTERED2ROPTIMIZEDMATRIX_H
#define FQHEMPSCLUSTERED2ROPTIMIZEDMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class RealMatrix;
class RealSymmetricMatrix;
class BosonOnDiskShort;
class AbstractArchitecture;

class FQHEMPSClustered2ROptimizedMatrix : public FQHEMPSClustered2RMatrix
{

  friend class FQHEMPSEvaluateCFTOperation;

 public:
  
  // default constructor 
  //
  FQHEMPSClustered2ROptimizedMatrix();

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool useRational = true, 
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool useRational = true,
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2ROptimizedMatrix(int pLevel, int nbrBMatrices, char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, 
			   double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSClustered2ROptimizedMatrix();
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName ();

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture);

  // get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
  //
  // cftSector = index of the CFT sector
  // return value = Q sector shift
  virtual int GetQValueCFTSectorShift(int cftSector);

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

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // compute the charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ);

  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

 protected:

  // compute the linearized index of a block of the matrices with fixed charge and p-level values
  //
  // chargedPartitionIndex = index of the partition in the charge sector
  // nbrCharges = total number of partitions at the current level
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
  // globalIndexShift = index of the first state at the considered level 
  // return value = linearized index
  virtual int Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
				      int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);


};

  
// get the name describing the B matrices 
// 
// return value = name 

inline char* FQHEMPSClustered2ROptimizedMatrix::GetName ()
{
  return this->BMatrixOutputName;
}

// compute the linearized index of a block of the matrices with fixed charge and p-level values
//
// chargedPartitionIndex = index of the partition in the charge sector
// nbrCharges = total number of partitions at the current level
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// globalIndexShift = index of the first state at the considered level 
// return value = linearized index

inline int FQHEMPSClustered2ROptimizedMatrix::Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
								      int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift)
{
  return (globalIndexShift + (chargedPartitionIndex + chargeSectorDimension * ((fieldIndex * nbrIdentityDescendant) + neutralPartitionIndex)));
}



#endif
