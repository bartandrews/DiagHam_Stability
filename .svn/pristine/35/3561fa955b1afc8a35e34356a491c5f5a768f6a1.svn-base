////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of MPS matrix for the PH-Pfaffian state               //
//                                                                            //
//                        last modification : 19/05/2017                      //
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


#ifndef FQHEMPSPHPFAFFIANMATRIX_H
#define FQHEMPSPHPFAFFIANMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class RealMatrix;
class RealSymmetricMatrix;
class BosonOnDiskShort;
class AbstractArchitecture;

class FQHEMPSPHPfaffianMatrix : public FQHEMPSClustered2RMatrix
{

 protected:

  friend class FQHEMPSEvaluateCFTOperation;

  // shift to apply to the p level to get a positive index
  int PLevelShift;


 public:
  
  // default constructor 
  //
  FQHEMPSPHPfaffianMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, const char* cftDirectory, bool bosonicVersion = false, bool useRational = true,
			  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0,
			  AbstractArchitecture* architecture = 0);


  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, const char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0);

  // destructor
  //
  ~FQHEMPSPHPfaffianMatrix();
  
  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (const char* cftDirectory, AbstractArchitecture* architecture);

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

  // get the array where the site-dependent matrices for the geometry are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseRealMatrix** GetSphereSiteDependentMatrices(int nbrFluxQuanta);

  // compute the CFT sector, the level and the charge index of a given matrix index
  //
  // index = matrix index
  // cftSector = reference on the CFT sector
  // pLevel = reference on the level
  // qValue = reference on the charge index
  // chargeSectorLevel = reference on the U(1) charge sector level
  // neutralSectorLevel = reference on the neutral sector level
  // chargeSectorIndex = index of the U(1) charge sector within its current level
  // neutralSectorIndex = index of the neutral sector within its current level
  virtual void GetCFTSectorChargeAndPLevelFromMatrixIndex(int index, int& cftSector, int& pLevel, int& qValue, int& chargeSectorLevel, 
							  int& neutralSectorLevel, int& chargeSectorIndex, int& neutralSectorIndex);

 protected:

  // compute the various arrays required to convert from quantum numbers and local indices to a global linearized index
  //
  // return value = dimension of the B matrix
  virtual int ComputeLinearizedIndexArrays();

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // pLevel = current total level
  // cftSector = CFT sector
  // qValue = charge index (i.e. Q)
  // chargeSectorLevel = level for the charge sector
  // chargeSectorIndex = index of the charge sector
  // cftSectorIndex = index within the current CFT sector and corresponding level (pLevel - chargeSectorLevel)
  // return value = linearized index
  virtual int Get2RMatrixIndexV2(int pLevel, int cftSector, int qValue, 
				 int chargeSectorLevel, int chargeSectorIndex, int cftSectorIndex);


};

  
// compute the linearized index of the B matrix for the (k=2,r) clustered states
//
// pLevel = current total level
// cftSector = CFT sector
// qValue = charge index (i.e. Q)
// chargeSectorLevel = level for the charge sector
// chargeSectorIndex = index of the charge sector
// cftSectorIndex = index within the current CFT sector and corresponding level (pLevel - chargeSectorLevel)
// return value = linearized index

inline int FQHEMPSPHPfaffianMatrix::Get2RMatrixIndexV2(int pLevel, int cftSector, int qValue, 
						       int chargeSectorLevel, int chargeSectorIndex, int cftSectorIndex)
{
  return (this->StartingIndexPerPLevelCFTSectorQValueU1Sector[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][chargeSectorLevel] 
	  + (this->NeutralSectorDimension[cftSector][this->PLevelShift + chargeSectorLevel - pLevel] * chargeSectorIndex) + cftSectorIndex);
}


// compute the CFT sector, the level and the charge index of a given matrix index
//
// index = matrix index
// sector = reference on the CFT sector
// pLevel = reference on the level
// qValue = reference on the charge index
// chargeSectorLevel = reference on the U(1) charge sector level
// neutralSectorLevel = reference on the neutral sector level
// chargeSectorIndex = index of the U(1) charge sector within its current level
// neutralSectorIndex = index of the neutral sector within its current level

inline void FQHEMPSPHPfaffianMatrix::GetCFTSectorChargeAndPLevelFromMatrixIndex(int index, int& cftSector, int& pLevel, int& qValue, int& chargeSectorLevel, 
										 int& neutralSectorLevel, int& chargeSectorIndex, int& neutralSectorIndex)
{
  pLevel = 0;
  cftSector = 0;
  qValue = 0;
  chargeSectorLevel = 0;
  neutralSectorLevel = 0;
  chargeSectorIndex = 0;
  neutralSectorIndex = 0;
  int MinIndexDistance = this->RealBMatrices->GetNbrRow();
  for (int p = 0; p <= this->PLevel; ++p)
    {
      for (int x = 0; x < this->NbrCFTSectors; ++x)
	{
	  for (int q = this->NInitialValuePerPLevelCFTSector[p][x]; q <= this->NLastValuePerPLevelCFTSector[p][x]; ++q)
	    {
	      int TmpIndexDistance = index - this->StartingIndexPerPLevelCFTSectorQValue[p][x][q - this->NInitialValuePerPLevelCFTSector[p][x]];
	      if ((TmpIndexDistance >= 0) && (TmpIndexDistance <= MinIndexDistance))
		{
		  MinIndexDistance = TmpIndexDistance;
		  pLevel = p;
		  cftSector = x;
		  qValue = q;
		}
	    }
	}     
    }  	
  MinIndexDistance = this->RealBMatrices->GetNbrRow();
  int TmpMinChargePSector = pLevel - this->PLevelShift;
  if (TmpMinChargePSector < 0)
    {
      TmpMinChargePSector = 0;
    }
  int TmpMaxChargePSector = pLevel;
  if (TmpMaxChargePSector > this->PLevelShift)
    {
      TmpMaxChargePSector = this->PLevelShift;
    }
  for (int TmpChargeLevel = TmpMinChargePSector; TmpChargeLevel <= TmpMaxChargePSector; ++TmpChargeLevel)
    {
      int TmpIndexDistance = index - this->StartingIndexPerPLevelCFTSectorQValueU1Sector[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][TmpChargeLevel];
      if ((TmpIndexDistance >= 0) && (TmpIndexDistance < MinIndexDistance))
	{
	  if (this->NeutralSectorDimension[cftSector][this->PLevelShift + TmpChargeLevel - pLevel] > 0)
	    {
	      MinIndexDistance = TmpIndexDistance;
	      chargeSectorLevel = TmpChargeLevel;
	    }
	}
    }
  index -= this->StartingIndexPerPLevelCFTSectorQValueU1Sector[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][chargeSectorLevel];
  neutralSectorLevel = (this->PLevelShift + chargeSectorLevel) - pLevel;
  chargeSectorIndex = index / this->NeutralSectorDimension[cftSector][neutralSectorLevel];  
  neutralSectorIndex = index % this->NeutralSectorDimension[cftSector][neutralSectorLevel];		  		    
}


#endif
