////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract MPS matrix for the FQHE                 //
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


#ifndef ABSTRACTFQHEMPSMATRIX_H
#define ABSTRACTFQHEMPSMATRIX_H


#include "config.h"
#include "MathTools/Complex.h" 
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"

#include <fstream>


using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

class AbstractFQHEMPSMatrix
{

 protected:

  // number of B matrices
  int NbrBMatrices;

  // arrar that describes all the physical indices
  unsigned long* PhysicalIndices;

  // array where the B matrices are stored (for real matrices)
  SparseRealMatrix* RealBMatrices;

  // array where the B matrices are stored (for complex matrices)
  SparseComplexMatrix* ComplexBMatrices;

  // array where the B matrices for quasiholes are stored
  SparseRealMatrix* QuasiholeBMatrices;

  // true if B matrix have to be evaluated on the torus geometry
  bool TorusFlag;

 public:
  
  // default constructor 
  //
  AbstractFQHEMPSMatrix();

  // destructor
  //
  ~AbstractFQHEMPSMatrix();
  
  // save the matrices 
  // 
  // fileName = name of the file where the matrices have to be stored
  // return value = true if no error occurred  
  virtual bool SaveMatrices (const char* fileName);

  // load the matrices 
  // 
  // fileName = name of the file where the matrices are stored
  // return value = true if no error occurred  
  virtual bool LoadMatrices (const char* fileName);

  // get the number of B matrices
  //
  // return value = number of B matrices
  virtual int GetNbrMatrices();

  // get the number of orbitals that associated to a set of B matrices
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // get the MPO bond dimension
  //
  // return value = MPO bond dimension
  virtual int GetBondDimension();

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseRealMatrix* GetMatrices();

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseComplexMatrix* GetComplexMatrices();

  // get the array where the site-dependent matrices are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseRealMatrix** GetSiteDependentMatrices(int nbrFluxQuanta);

  // get the array where the site-dependent matrices are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseComplexMatrix** GetSiteDependentComplexMatrices(int nbrFluxQuanta);

  // get the array where the site-dependent matrices for the geometry are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseRealMatrix** GetSphereSiteDependentMatrices(int nbrFluxQuanta);

  // get the array where the site-dependent matrices for the geometry are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseComplexMatrix** GetSphereSiteDependentComplexMatrices(int nbrFluxQuanta);

  // get the edge matrix for localized quasiholes, with normal ordering
  //
  // nbrQuasiholes = number of quasiholes
  // quasiholePositions = quasihole positions (for cylinder, positions have to be expressed in perimeter units)
  // return value = pointer to the edge matrix
  virtual SparseComplexMatrix* GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions);
  
  // test if the MPS is written for the torus geometry
  //
  // return value = true if the MPS is written for the torus geometry
  virtual bool IsTorus();

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
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the (k,r) exclude principle satisfied by the root configuration
  // 
  // pauliK = maximum number of particles in the pauliR consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual void GetKRExclusionPrinciple(int& pauliK, int& pauliR); 

  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

  // get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
  // 
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // statistics = true if we are dealing with fermions
  // return value = minimum ky momentum 
  virtual int GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics);
  
  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // get the MPS truncation level
  //
  // return value = truncation level
  virtual int GetTruncationLevel();

  // get the number of CFT sectors invloved on the MPS
  //
  // return value = number of CFT sectors
  virtual int GetNbrCFTSectors();

  // get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
  //
  // cftSector = index of the CFT sector
  // return value = Q sector shift
  virtual int GetQValueCFTSectorShift(int cftSector);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual RealMatrix ExtractBlock(RealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // cftSector1 = CFT sector of the blck left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // cftSector2 = CFT sector of the blck right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int cftSector1, int q1, int pLevel2, int cftSector2, int q2);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // cftSector1 = CFT sector of the blck left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // cftSector2 = CFT sector of the blck right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual RealMatrix ExtractBlock(RealMatrix& matrix, int pLevel1, int cftSector1, int q1, int pLevel2, int cftSector2, int q2);

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

  // get the charge index range at a given truncation level
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

  // compute the global charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index

  virtual void ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ);

  // compute P, N from the linearized index of the B matrix for the Laughlin states
  //
  // index = linearized index
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  virtual void GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex);

  // compute the level and the charge index of a given matrix index
  //
  // index = matrix index
  // pLevel = reference on the level
  // qValue = reference on the charge index
  virtual void GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue);

  // compute the CFT sector, the level and the charge index of a given matrix index
  //
  // index = matrix index
  // sector = reference on the CFT sector
  // pLevel = reference on the level
  // qValue = reference on the charge index
  virtual void GetCFTSectorChargeAndPLevelFromMatrixIndex(int index, int& sector, int& pLevel, int& qValue);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the extra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // get the parent MPS matrices if the current MPS matrices have ones
  //
  // return value = pointer to the parent MPS matrices
  virtual AbstractFQHEMPSMatrix* GetParentMPSMatrices();

  // get the array that gives the index of each entry of the current B matrix within its parent B matrix
  //
  // return value = index mapping array
  virtual int* GetIndexMappingArray();

  // get the array of physical indices
  //
  // return value  = array of physical indices
  virtual unsigned long* GetPhysicalIndices();

  // get a given physical index
  //
  // index = index to retrieve
  //  configuration = array where the description of the physical index will be stored
  virtual void GetPhysicalIndex(int index, unsigned long* configuration);

  // print a given physical index
  //
  // str = reference on the output stream
  // index = integer associated to the  physical index 
  // return value = reference on the output stream
  virtual ostream& PrintPhysicalIndex(ostream& str, int index);

  // get a given physical index at a given orbital for the site-dependent MPS
  //
  // orbitalIndex = orbital index 
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  virtual void GetSiteDependentPhysicalIndex(int orbitalIndex, int index, unsigned long* configuration);

  // print a given state of the auxiliary space
  //
  // str = reference on the output stream
  // index = index of the state
  // return value = reference on the output stream
  virtual ostream& PrintAuxiliarySpaceState(ostream& str, int index);
 
  // get the label associated to a given state of the auxiliary space
  //
  // index = auxiliary space index
  // return value = string containing the label
  virtual char* GetAuxiliarySpaceLabel(int index);

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

  // load the specific informations from the file header
  // 
  // file = reference on the input file stream
  // return value = true if no error occurred  
  virtual bool LoadHeader (ifstream& file);

  // save the specific informations to the file header 
  // 
  // file = reference on the output file stream
  // return value = true if no error occurred  
  virtual bool SaveHeader (ofstream& file);

};

// get the number of B matrices
//
// return value = number of B matrices

inline int AbstractFQHEMPSMatrix::GetNbrMatrices()
{
  return this->NbrBMatrices;
}

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseRealMatrix* AbstractFQHEMPSMatrix::GetMatrices()
{
  return this->RealBMatrices;
}

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseComplexMatrix* AbstractFQHEMPSMatrix::GetComplexMatrices()
{
  return this->ComplexBMatrices;
}

// get the number of orbitals that associated to a set of B matrices
//
// return value = number oforbitals

inline int AbstractFQHEMPSMatrix::GetNbrOrbitals()
{
  return 1;
}

// get the maximum occupation per orbital
//
// return value = aximum occupation per orbital

inline int AbstractFQHEMPSMatrix::GetMaximumOccupation()
{
  return 1;
}

// compute P, N from the linearized index of the B matrix for the Laughlin states
//
// index = linearized index
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector

inline void AbstractFQHEMPSMatrix::GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex)
{
  cout << "Dummy. " << endl;
}

// compute the level and the charge index of a given matrix index
//
// index = matrix index
// pLevel = reference on the level
// qValue = reference on the charge index

inline void AbstractFQHEMPSMatrix::GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue)
{
  pLevel = -1;
  qValue = -1;
}


// compute the CFT sector, the level and the charge index of a given matrix index
//
// index = matrix index
// sector = reference on the CFT sector
// pLevel = reference on the level
// qValue = reference on the charge index

inline void AbstractFQHEMPSMatrix::GetCFTSectorChargeAndPLevelFromMatrixIndex(int index, int& sector, int& pLevel, int& qValue)
{
  sector = 0;
  this->GetChargeAndPLevelFromMatrixIndex(index, pLevel, qValue);
}


// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int AbstractFQHEMPSMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return 1;
}

// get the MPS truncation level
//
// return value = truncation level

inline int AbstractFQHEMPSMatrix::GetTruncationLevel()
{
  return 0;
}

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int AbstractFQHEMPSMatrix::GetNbrCFTSectors()
{
  return 1;
}

// get the array of physical indices
//
// return value  = array of physical indices

inline unsigned long* AbstractFQHEMPSMatrix::GetPhysicalIndices()
{
  return this->PhysicalIndices;
}

// test if the MPS is written for the torus geometry
//
// return value = true if the MPS is written for the torus geometry

inline bool AbstractFQHEMPSMatrix::IsTorus()
{
  return this->TorusFlag;
}

// get the MPO bond dimension
//
// return value = MPO bond dimension

inline int AbstractFQHEMPSMatrix::GetBondDimension()
{
  if (this->RealBMatrices != 0)
    {
      return this->RealBMatrices[0].GetNbrRow();
    }
  else
    {
      return this->ComplexBMatrices[0].GetNbrRow();
    }
}

#endif
