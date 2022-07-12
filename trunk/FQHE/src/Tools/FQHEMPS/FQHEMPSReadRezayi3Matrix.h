////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
//                                                                            //
//                        last modification : 30/10/2012                      //
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


#ifndef FQHEMPSREADREZAYI3MATRIX_H
#define FQHEMPSREADREZAYI3MATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"


class FQHEMPSReadRezayi3Matrix : public FQHEMPSClustered2RMatrix
{

 protected:


 public:
  
  // default constructor 
  //
  FQHEMPSReadRezayi3Matrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
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
  FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool bosonicVersion = false, bool useRational = true, 
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			   bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0,
			   AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
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
  FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion = false, bool useRational = true, 
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
			   bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0, AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0, 
			   bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0);

  // destructor
  //
  ~FQHEMPSReadRezayi3Matrix();
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName ();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory = 0, AbstractArchitecture* architecture = 0);

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

  // get the matrix that into account the Jordan Wigner string on the torus geometry
  //
  // nbrFermions = number of fermions in the system
  // return value = corresponding matrix
  virtual SparseRealMatrix GetTorusStringMatrix(int nbrFermions);

  // get the auxiliary space indices that are related to a given topological scetor
  //
  // topologicalSector = index of the topological sector to select
  // nbrIndices = reference on the integer that will be set to the number of indices
  // return value = array that contains the auxiliary space indices related to the selected topological sector
  virtual int* GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices);
  
 protected:

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  // nbrCharges = total number of charge indices
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi_{+1}, 3 for psi_{-1}, 5 for W_-3)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendants of the identity at the current level
  // nbrPsiDescendant = number of linearly independent descendants of the psi_{+1} at the current level
  // globalIndexShift = index of the first state at the considered level
  // return value = linearized index
  int GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, 
				 int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int nbrPsiDescendant, 
				 int globalIndexShift);


};


// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSReadRezayi3Matrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return 5;
}

// compute the linearized index of the B matrix for the (k=2,r) clustered states
//
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector
// nbrCharges = total number of charge indices
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// nbrPsiDescendant = number of linearly independent descendants of the psi_{+1} at the current level
// globalIndexShift = index of the first state at the considered level
// return value = linearized index

inline int FQHEMPSReadRezayi3Matrix::GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, 
								int nbrCharges, int chargeSectorDimension, 
								int fieldIndex, int neutralPartitionIndex, 
								int nbrIdentityDescendant, int nbrPsiDescendant, 
								int globalIndexShift)
{
  return ((((((nbrIdentityDescendant * (fieldIndex & 1)) + (nbrPsiDescendant * (fieldIndex >> 1))) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}

#endif
